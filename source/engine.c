#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi (3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348)
#define two_pi (6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696)
#define polyphony (100)
#define max_harms (30)
#define maximum_decay (5)
#define noise_ncomp (20)

float dynamics_factor = 0.5f, noise_value = 100.f, K_multiplier = 1.0f, delay_note = 0.02f, noise_decay = 1.0f, Q_multiplier = 1.0f, power_noise = 0.3f, stereo = 1.0f;
double t_current = 0.;

int count= 0;


double Mod(double x, double y){
	if (y == 0.)
		return x;

	double m= x - y * floor(x/y);
	if (y > 0){
		if (m >= y)
			return 0;

		if (m<0 ){
			if (y+m == y)
				return 0;
			else
				return y+m;
		}
	}
	else{
		if (m <= y)
			return 0;

		if (m > 0){
			if (y+m == y)
				return 0;
			else
				return y+m;
		}
	}

	return m;
}


void set_parameters(float dynamics, float noise_level, float volume, float K, float delay, float noise_t, float Q_mult, float noise_exp){
	dynamics_factor = dynamics;
	noise_value     = noise_level;
	K_multiplier    = K;
	delay_note      = delay;
	noise_decay     = noise_t;
	Q_multiplier    = Q_mult;
	power_noise     = noise_exp;
}


static inline double random(){
	return (double)rand()/(double)(RAND_MAX);
}

static inline float get_oscillator_value(float t, float t_init, float A, float f, float t0, float delay, float phase){

	// else if(t >= t0*4 + t_init)
	// 	return 0.0f;

	if(t0 == 0.0f)
		return 0.0f;
	else{
		if(t < t_init + delay){
			return 0.0f;
			// printf("juajuaja %f\n", t - t_init - delay);
		}

		else
			return A*sinf(2*pi*f*(t - t_init - delay) + phase)*expf(-(t - t_init - delay)/t0);
	}
}


struct harmonic_info{
	int min_key_value;
	int key_amount;

	int   *nharms;
	float *freqs, *value_l, *value_r, *f1, *k, *b1, *b2;
	float *t_decays, *delays;
	float *amp_l, *amp_r;
};

struct harmonic_info harmonic_content;


struct guitar_sound{
	// status == 0 => inactive
	// status == 1 => active
	// status == 2 => quick fade out after key release

	int status, nharms, priority, key_id, key_idx;
	float t_init, time_keyoff, total_decay, keyoff_decay, delay, b1, b2, k, f1, *amp_l, *amp_r, *t_decays, *freqs, *phase, *delays;
};


struct strum_noise{
	int status, priority;
	float *freqs, *Qs, *t_decay1, *t_decay2, total_decay, t_init, low_f1, low_f2, amp1, amp2, *amplitude, *x2, *x1, *y1, *y2, *B0, *B1, *B2, *A1, *A2;
};

struct guitar_sound oscillators[polyphony];
struct strum_noise  noise_osc[polyphony];

float *harmonic_multiplier;


void disable_oscillators(){
	for(unsigned int i = 0; i < polyphony; i++){
		oscillators[i].status = 0;
		noise_osc[i].status = 0;
	}
}

void reset_timer(){
	t_current = 0.f;
	disable_oscillators();
}

static inline float model_k(int n, float f0, float K){
	// return f0*sqrtf(1.0f + K*K*n*n)*n;
	return f0*sqrtf(1.0f + K*K*n*n)*n/sqrtf(1.0f + K*K);
}

static inline float model_t0(int n, float b1, float b2){
	return 1.0f/(b1 + n*n*b2);

}
void initialize_oscillators(){
	for(unsigned int i = 0; i < polyphony; i++){
		oscillators[i].amp_l    = (float *) malloc(max_harms*sizeof(float));
		oscillators[i].amp_r    = (float *) malloc(max_harms*sizeof(float));
		oscillators[i].freqs    = (float *) malloc(max_harms*sizeof(float));
		oscillators[i].t_decays = (float *) malloc(max_harms*sizeof(float));
		oscillators[i].phase    = (float *) malloc(2*max_harms*sizeof(float));
		oscillators[i].delays   = (float *) malloc(2*max_harms*sizeof(float));

		noise_osc[i].freqs     = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].Qs        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].t_decay1  = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].t_decay2  = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].amplitude = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].B0        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].B1        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].B2        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].A1        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].A2        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].x1        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].x2        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].y1        = (float *) malloc(noise_ncomp*sizeof(float));
		noise_osc[i].y2        = (float *) malloc(noise_ncomp*sizeof(float));


	}
}


void load_harmonic_multipliers(float *multipliers, unsigned int nharms){
	for(unsigned int i = 0; i < nharms; i++)
		harmonic_multiplier[i] = multipliers[i];
}


void initialize_harmonic_info(int min_key_value, unsigned int key_amount, float *freqs, float *t_decays, float *amplitudes_left, float *amplitudes_right, float *b1s, float *b2s, float *ks, float *f1s, float *values_l, float *values_r, int *nharms){
	float period;
	harmonic_content.min_key_value = min_key_value;
	harmonic_content.key_amount    = key_amount;

	harmonic_multiplier         = (float *) malloc(max_harms*sizeof(float));
	for(unsigned int j = 0; j < max_harms; j++) harmonic_multiplier[j] = 1.0f;

	harmonic_content.freqs      = (float *) malloc(key_amount*max_harms*sizeof(float));
	harmonic_content.t_decays   = (float *) malloc(key_amount*max_harms*sizeof(float));
	harmonic_content.amp_l      = (float *) malloc(key_amount*max_harms*sizeof(float));
	harmonic_content.amp_r      = (float *) malloc(key_amount*max_harms*sizeof(float));
	harmonic_content.delays     = (float *) malloc(2*key_amount*max_harms*sizeof(int));
	harmonic_content.b1         = (float *) malloc(key_amount*sizeof(float));
	harmonic_content.b2         = (float *) malloc(key_amount*sizeof(float));
	harmonic_content.k          = (float *) malloc(key_amount*sizeof(float));
	harmonic_content.f1         = (float *) malloc(key_amount*sizeof(float));
	harmonic_content.value_l    = (float *) malloc(key_amount*sizeof(float));
	harmonic_content.value_r    = (float *) malloc(key_amount*sizeof(float));
	harmonic_content.nharms     = (int *)   malloc(key_amount*sizeof(int));

	for(unsigned int i = 0; i < key_amount; i++){
		printf("to load key: %d\n", i);
		harmonic_content.b1[i]      = b1s[i];
		harmonic_content.b2[i]      = b2s[i];
		harmonic_content.k[i]       = ks[i];
		harmonic_content.f1[i]      = f1s[i];
		harmonic_content.value_l[i] = values_l[i];
		harmonic_content.value_r[i] = values_r[i];
		harmonic_content.nharms[i]  = nharms[i];

		for(unsigned int j = 0; j < max_harms; j++){
			harmonic_content.freqs[i*max_harms + j]    = freqs[i*max_harms + j];
			harmonic_content.t_decays[i*max_harms + j] = t_decays[i*max_harms + j];
			harmonic_content.amp_l[i*max_harms + j]    = amplitudes_left[i*max_harms + j];
			harmonic_content.amp_r[i*max_harms + j]    = amplitudes_right[i*max_harms + j];

			period = 1.0f/harmonic_content.freqs[i*max_harms + j];
			harmonic_content.delays[2*(i*max_harms + j) + 0] = period*random()*stereo;
			harmonic_content.delays[2*(i*max_harms + j) + 1] = period*random()*stereo;
		}
	}
	printf("loaded 'successfully'\n");
}


float noise_envelope(float t, float delta1, float delta2, float t_centre){
	if(t <= t_centre - delta1)
		return 0.f;
	else if(t <= t_centre)
		return powf((t - t_centre + delta1)/delta1, 2);
	else if(t <= t_centre + delta2)
		return expf(-(t - t_centre)/delta2);
	else
		return 0.f;
}

void initialize_strum_noise(int osc_id, int key_id, float time, float velocity, unsigned int samplerrate){
	float rand, f_centre, K, norm, Q, b0, b1, b2, a1, a2;
	// noise_osc[osc_id].t_center    = 0.02f;
	noise_osc[osc_id].priority    = polyphony;
	noise_osc[osc_id].status      = 1;
	noise_osc[osc_id].t_init      = time;
	noise_osc[osc_id].total_decay = 0.1f;

	float factor = (key_id - 18)/(42.0f);

	noise_osc[osc_id].low_f1 = 83.0f  + random()*5.0f;
	noise_osc[osc_id].low_f2 = 170.0f + random()*10.0f;
	noise_osc[osc_id].amp1   = sqrtf(0.3f)*factor;
	noise_osc[osc_id].amp2   = sqrtf(0.2f)*factor;


	// float max_freq = logf(12000.0f/2000)/logf(10.f)*factor + logf(2000)/logf(10.f);
	float max_freq = logf(12000.0f/2000)/logf(10.f)*factor + logf(2000)/logf(10.f);
//	float min_freq = logf(1000.0f/300)/logf(10.f)*factor + logf(300)/logf(10.f);
	// float f_centres[10] = {250., 300.f, 500.f, 1500.f, 1800.f, 6000.f, 7000.f, 8000.f, 9000.f, 12000.f};
	float f_centres[20] = {250., 300.f, 500.f, 1500.f, 1800.f, 6000.f, 7000.f, 8000.f, 9000.f, 12000.f, 600.f, 650.f, 900.f, 1200.f, 2000.f, 2500.f, 2800.f, 3500.f, 4500.f, 5500.f};
	float Qs[10] = {50.f, 50.f, 30.f, 50.f, 50.f, 40.f, 30.f, 10.f, 5.f, 5.f};
	float min_freq = logf(300)/logf(10.f);

	for(unsigned int i = 0; i < noise_ncomp; i++){
		rand = random();

		f_centre = f_centres[i]*(9.5 + random())/10.0;
		Q        = 50*Q_multiplier;


		noise_osc[osc_id].amplitude[i] = 6*powf(400.0f/f_centre, 2 - 1.2)*(factor + .1f)/2.*(velocity*velocity)*rand;

		noise_osc[osc_id].t_decay1[i]  = 0.012f*powf(1000.0f/f_centre, 0.7f)*noise_decay;
		noise_osc[osc_id].t_decay2[i]  = 0.008f*powf(1000.0f/f_centre, 0.7f)*noise_decay;

		noise_osc[osc_id].Qs[i]        = Q;
		noise_osc[osc_id].freqs[i]     = f_centre;

		K = tanf(pi*f_centre/samplerrate);

		norm = 1.0f / (1.0f + K / Q + K * K);
		b0   = K / Q * norm;
		b1   = 0.0f;
		b2   = -b0;
		a1   = 2 * (K * K - 1) * norm;
		a2   = (1 - K / Q + K * K) * norm;

		noise_osc[osc_id].B0[i] = b0;
		noise_osc[osc_id].B1[i] = b1;
		noise_osc[osc_id].B2[i] = b2;
		noise_osc[osc_id].A1[i] = a1;
		noise_osc[osc_id].A2[i] = a2;
		noise_osc[osc_id].x2[i] = 0.0f;
		noise_osc[osc_id].x1[i] = 0.0f;
		noise_osc[osc_id].y1[i] = 0.0f;
		noise_osc[osc_id].y2[i] = 0.0f;
	}
}




float strum_noise_generator(float *result, unsigned int N, unsigned int samplerrate){
	float t, dt = 1.0f/samplerrate, x1, x2, y1, y2, b0, b1, b2, a1, a2, x0, y0, t_centre = 0.012f, A, env, A1, A2, t_decay1, t_decay2;
	for(unsigned int osc_id = 0; osc_id < polyphony; osc_id++){
		if(noise_osc[osc_id].status == 1){
			A1 = noise_osc[osc_id].amp1*noise_value*0;
			A2 = noise_osc[osc_id].amp2*noise_value*0;
			for(unsigned int i = 0; i < noise_ncomp; i++){
				x1 = noise_osc[osc_id].x1[i];
				x2 = noise_osc[osc_id].x2[i];
				y1 = noise_osc[osc_id].y1[i];
				y2 = noise_osc[osc_id].y2[i];
				b0 = noise_osc[osc_id].B0[i];
				b1 = noise_osc[osc_id].B1[i];
				b2 = noise_osc[osc_id].B2[i];
				a1 = noise_osc[osc_id].A1[i];
				a2 = noise_osc[osc_id].A2[i];

				t_decay1 = noise_osc[osc_id].t_decay1[i];
				t_decay2 = noise_osc[osc_id].t_decay2[i];
				A        = noise_osc[osc_id].amplitude[i]*noise_value;

				for(unsigned int k = 0; k < N; k++){
					t = t_current + k*dt;

					x0 = random();
					y0 = b0*x0 + b1*x1 + b2*x2 - a1*y1 - a2*y2;

					env  = noise_envelope(t, t_decay1, t_decay2, t_centre + noise_osc[osc_id].t_init);
					result[2*k + 0] += y0*A*env;
					result[2*k + 1] += y0*A*env;

					x2 = x1;
					x1 = x0;
					y2 = y1;
					y1 = y0;
				}
				noise_osc[osc_id].x2[i] = x2;
				noise_osc[osc_id].x1[i] = x1;
				noise_osc[osc_id].y2[i] = y2;
				noise_osc[osc_id].y1[i] = y1;


			}
			for(unsigned int k = 0; k < N; k++){
				t = t_current + k*dt;
				result[2*k + 0] += get_oscillator_value(t, noise_osc[osc_id].t_init, A1, noise_osc[osc_id].low_f1, 0.1f,  t_centre, 0.);
				result[2*k + 1] += get_oscillator_value(t, noise_osc[osc_id].t_init, A1, noise_osc[osc_id].low_f1, 0.1f,  t_centre, 0.);
				result[2*k + 0] += get_oscillator_value(t, noise_osc[osc_id].t_init, A2, noise_osc[osc_id].low_f2, 0.01f, t_centre, 0.);
				result[2*k + 1] += get_oscillator_value(t, noise_osc[osc_id].t_init, A2, noise_osc[osc_id].low_f2, 0.01f, t_centre, 0.);
			}
		}

	}



}


void initialize_guitar_sound(int osc_id, int key_id, float time, float velocity){
	int key_idx = key_id - harmonic_content.min_key_value;
	int nharms = harmonic_content.nharms[key_idx];
	float f0_theoretical = 92.49860567790869*powf(1.0594630943592953, key_id - 18), dynamics, period;

	oscillators[osc_id].status   = 1;
	oscillators[osc_id].priority = polyphony;
	oscillators[osc_id].key_id   = key_id;
	oscillators[osc_id].key_idx  = key_idx;
	oscillators[osc_id].t_init   = time;

	oscillators[osc_id].b1 = harmonic_content.b1[key_idx];
	oscillators[osc_id].b2 = harmonic_content.b2[key_idx];
	oscillators[osc_id].k  = harmonic_content.k[key_idx];
	oscillators[osc_id].f1 = harmonic_content.f1[key_idx]*0 + f0_theoretical;
	oscillators[osc_id].nharms = harmonic_content.nharms[key_idx];
	oscillators[osc_id].delay  = delay_note;

	for(unsigned int i = 0; i < nharms; i++){
		dynamics = fmaxf(1 + (velocity - 0.3f)*i*dynamics_factor/nharms, 0.f);
		// oscillators[osc_id].freqs[i]    = harmonic_content.freqs[key_idx*max_harms + i];
		oscillators[osc_id].freqs[i]    = model_k(i + 1, f0_theoretical, oscillators[osc_id].k*K_multiplier);
		oscillators[osc_id].t_decays[i] = harmonic_content.t_decays[key_idx*max_harms + i];
		oscillators[osc_id].amp_l[i]    = velocity*harmonic_content.amp_l[key_idx*max_harms + i]*dynamics;
		oscillators[osc_id].amp_r[i]    = velocity*harmonic_content.amp_r[key_idx*max_harms + i]*dynamics;

		period                          = 1.0f/oscillators[osc_id].freqs[i];

		oscillators[osc_id].delays[2*i + 0]   = delay_note + harmonic_content.delays[2*(key_idx*max_harms + i) + 0];
		oscillators[osc_id].delays[2*i + 1]   = delay_note + harmonic_content.delays[2*(key_idx*max_harms + i) + 1];

		oscillators[osc_id].phase[2*i + 0]    = oscillators[osc_id].delays[2*i + 0]/period*0;
		oscillators[osc_id].phase[2*i + 1]    = oscillators[osc_id].delays[2*i + 1]/period*0;
	}
	oscillators[osc_id].total_decay = 1.0f/oscillators[osc_id].b1*maximum_decay;
}

int check_lowest_priority(){
	int minimum_priority = polyphony, osc_id = 0;

	for(unsigned int i = 0; i < polyphony; i++){
		if(oscillators[i].priority < minimum_priority){
			minimum_priority = oscillators[i].priority;
			osc_id = i;
		}
	}
	return osc_id;
}

int lowest_priority_noise(){
	int minimum_priority = polyphony, osc_id = 0;

	for(unsigned int i = 0; i < polyphony; i++){
		if(noise_osc[i].priority < minimum_priority){
			minimum_priority = noise_osc[i].priority;
			osc_id = i;
		}
	}
	return osc_id;
}
int get_oscillator_id_from_key(int key_id, int *found){
	int result = 0;
	*found = 0;
	for(unsigned int i = 0; i < polyphony; i++){
		if((oscillators[i].key_id == key_id)&(oscillators[i].status == 1)){
			result = i;
			*found = 1;
			break;
		}
	}
	return result;
}

void play_note(unsigned int key_id, float velocity, int samplerrate){
	int osc_id, found = 0;

	osc_id = get_oscillator_id_from_key(key_id, &found);
	if((found == 0)|(oscillators[osc_id].status != 2)){
		osc_id = check_lowest_priority();
	}
	initialize_guitar_sound(osc_id, key_id, t_current, velocity);

	osc_id = lowest_priority_noise();
	initialize_strum_noise(osc_id, key_id, t_current, velocity, samplerrate);

	for(unsigned int i = 0; i < polyphony; i++){
		oscillators[i].priority += -1;
		noise_osc[i].priority   += -1;
	}
}


void key_off(unsigned int key_id, unsigned int samplerrate){
	int osc_id, found;

	float dt = 1.0f/samplerrate, phase1, phase2, delay1, delay2;
	osc_id = get_oscillator_id_from_key(key_id, &found);
	if((found == 1)&(oscillators[osc_id].status == 1)){
		oscillators[osc_id].status = 2;
		oscillators[osc_id].priority    = polyphony;

		float time_ant = oscillators[osc_id].t_init;
		oscillators[osc_id].t_init = t_current;
		for(unsigned int i = 0; i < oscillators[osc_id].nharms; i++){
			delay1 = oscillators[osc_id].delays[2*i + 0];
			delay2 = oscillators[osc_id].delays[2*i + 1];

			phase1 = Mod(2*pi*oscillators[osc_id].freqs[i]*(t_current -  time_ant - delay1), two_pi);
			phase2 = Mod(2*pi*oscillators[osc_id].freqs[i]*(t_current - t_current), two_pi);
			oscillators[osc_id].phase[2*i + 0]    = phase1 - phase2;

			phase1 = Mod(2*pi*oscillators[osc_id].freqs[i]*(t_current -  time_ant - delay2), two_pi);
			phase2 = Mod(2*pi*oscillators[osc_id].freqs[i]*(t_current - t_current), two_pi);
			oscillators[osc_id].phase[2*i + 1]    = phase1 - phase2;


			oscillators[osc_id].amp_l[i]    *= expf(-(t_current - time_ant - delay1)/oscillators[osc_id].t_decays[i]);
			oscillators[osc_id].amp_r[i]    *= expf(-(t_current - time_ant - delay2)/oscillators[osc_id].t_decays[i]);
			oscillators[osc_id].t_decays[i]  = fminf(oscillators[osc_id].t_decays[i]/16.0f, 0.1f);

			oscillators[osc_id].delays[2*i + 0] = 0.f;
			oscillators[osc_id].delays[2*i + 1] = 0.f;
		}
		oscillators[osc_id].delay        = 0.f;
		oscillators[osc_id].total_decay *= 1.f/16;
	}
	osc_id = lowest_priority_noise();
	initialize_strum_noise(osc_id, key_id, t_current, 0.15f, samplerrate);
	
}

void check_status(float time){
	for(unsigned int i = 0; i < polyphony; i++){
		if((oscillators[i].status == 1)|(oscillators[i].status == 2)){
			if(time - oscillators[i].t_init >= oscillators[i].total_decay){
				oscillators[i].status = 0;
			}
		}
		if(noise_osc[i].status == 1){
			if(time - noise_osc[i].t_init >= noise_osc[i].total_decay){
				noise_osc[i].status = 0;
			}
		}
	}
}



void deallocate_memory(){
	for(unsigned int i = 0; i < polyphony; i++){
		free(oscillators[i].amp_l);
		free(oscillators[i].amp_r);
		free(oscillators[i].freqs);
		free(oscillators[i].t_decays);
	}
}



void fill_chunk(float *values, unsigned int N, unsigned int samplerrate){
	float time, dt = 1.0f/samplerrate, A_l, A_r, t0, f, f1, K, t_init, phase1, phase2, delay;
	int *status, j, nharms, active_oscillators = 0, sign;
	for(j = 0; j < N; j++){
		values[2*j + 0] = 0.f;
		values[2*j + 1] = 0.f;
	}
	// printf("Current time: %f\n", t_current);
	for(unsigned int osc_id = 0; osc_id < polyphony; osc_id++){
		if(oscillators[osc_id].status != 0){
			nharms = oscillators[osc_id].nharms;
			t_init = oscillators[osc_id].t_init;
			f1  = oscillators[osc_id].f1;
			K   = oscillators[osc_id].k*K_multiplier;
			delay = oscillators[osc_id].delay;
			// printf("Status: %d, K = %f\n", oscillators[osc_id].status, K);

			sign = -1;
			for(int k = 0; k < nharms; k++){
				A_l = oscillators[osc_id].amp_l[k]*harmonic_multiplier[k];
				A_r = oscillators[osc_id].amp_r[k]*harmonic_multiplier[k];
				t0  = oscillators[osc_id].t_decays[k];
				phase1 = oscillators[osc_id].phase[2*k + 0];
				phase2 = oscillators[osc_id].phase[2*k + 1];
				f   = oscillators[osc_id].freqs[k];
				// f   = model_k(k + 1, f1, K);
				sign *= (-1);
				for(j = 0; j < N; j++){
					time = t_current + j*dt;

					values[2*j + 0] += get_oscillator_value(time, t_init, A_l, f, t0, oscillators[osc_id].delays[2*k + 0], phase1);
					values[2*j + 1] += get_oscillator_value(time, t_init, A_r, f, t0, oscillators[osc_id].delays[2*k + 1], phase2);
				}
			}
			active_oscillators += 1;
		}
	}
	strum_noise_generator(values, N, samplerrate);

	printf("Active oscillators: %d  \r", active_oscillators);
	t_current += N*dt;

	check_status(t_current);
}

//gcc -c -Ofast -fPIC source/engine.c -o lib/engine.o -lm -fopenmp -std=c99
//gcc -shared lib/engine.o -o lib/engine.so -lm -fopenmp

