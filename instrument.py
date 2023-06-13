import os
import numpy as np
from scipy.io import wavfile
import engine_module as engine
import pyaudio


def load_info():

	name_keysavail = os.listdir('harms/')
	keys_avail = []
	for item in name_keysavail:
		item = item.replace("key_", "")
		item = item.replace(".npy", "")
		keys_avail.append(int(item))
	keys_avail = np.array(keys_avail)
	keys_avail.sort()

	minimum_key_id = keys_avail.min()

	amps_l, amps_r, freqs, t_decays, b1s, b2s, ks, f1s, values_l, values_r, nharms = [], [], [], [], [], [], [], [], [], [], []
	for key_id in keys_avail:
		print(key_id)
		info = np.load("harms/key_" + str(key_id) + ".npy", allow_pickle = True).item()

		amps_l.append(info["A"][0])
		amps_r.append(info["A"][1])
		freqs.append(info["f"])
		t0 = (info["t0"][0] + info["t0"][1])/2.
		t_decays.append(t0)
		b1s.append(info["b1"])
		b2s.append(info["b2"])
		ks.append(info["k"])
		f1s.append(info["f1"])
		values_l.append(info["volume"][0])
		values_r.append(info["volume"][1])
		nharms.append(info["nharms"])

	amps_l = np.concatenate(amps_l)
	amps_r = np.concatenate(amps_r)
	freqs  = np.concatenate(freqs)
	t_decays = np.concatenate(t_decays)
	b1s      = np.array(b1s)
	b2s      = np.array(b2s)
	f1s      = np.array(f1s)
	values_l = np.array(values_l)
	values_r = np.array(values_r)
	nharms   = np.array(nharms)

	engine.load_harmonics(minimum_key_id, len(keys_avail), freqs, t_decays, amps_l, amps_r, b1s, b2s, ks, f1s, values_l, values_r, nharms)


def set_multipliers(mults):
	multipliers = np.ones(30)
	for i in range(min(30, len(mults))):
		multipliers[i] = mults[i]

	engine.set_multipliers(multipliers)




load_info()

engine.initialize_engine()




if False:
	from pygame import midi

	midi.init()
	default_id = midi.get_default_input_id()
	midi_input = midi.Input(device_id=3)



class AudioFile:
	def __init__(self, params):
		""" Init audio stream """
		self.p = pyaudio.PyAudio()
		self.record = False
		self.params = params
		engine.set_parameters(self.params)

		self.stream = self.p.open(format=pyaudio.paFloat32,
					channels=2,
					rate=44100,
					output=True)

	def close(self):
		if self.record:
			self.audio = np.concatenate(self.audio)

			two_chanels = np.zeros([len(self.audio), 2])
			two_chanels[:, 0] = self.audio[:, 0]
			two_chanels[:, 1] = self.audio[:, 1]

			if not os.path.isdir("output"):
				os.system("mkdir output")

			wavfile.write("output/output.wav", data = two_chanels, rate = 44100)

			self.stream.close()
			self.p.terminate()

		if self.midi_record:
			data = {"idx": self.midi_time, "msg": self.midi_signal}
			np.save("midi_demo", data)


	def play(self, record = False, record_midi = False, chunk = 256, samplerrate = 44100, play_midi = False):
		self.record, self.midi_record  = record, record_midi
		if play_midi: self.midi_record = False

		self.audio, count, self.midi_signal, self.midi_time = [], 0, [], []

		if play_midi:
			data = np.load("midi_demo.npy", allow_pickle = True).item()
			midi_signal, midi_time, note_count = data["msg"], data["idx"], 0
		while True:
			if play_midi:
				if count == midi_time[note_count]:
					qwe = True
				else:
					qwe = False
			else:
				qwe = midi_input.poll()

			if qwe:
				if play_midi:
					mess = midi_signal[note_count]
					note_count += 1
				else:
					mess = midi_input.read(1)
					self.midi_signal.append(mess)
					self.midi_time.append(count)


				if mess[0][0][0] == 144:
					key_id = mess[0][0][1] - 42 + 18
					vel = mess[0][0][2]
					if (key_id >= 18) and ( key_id < 18 + 42):
						engine.play_note(key_id, vel/128, samplerrate)

					elif mess[0][0][1] == 35:
						self.close()
						break

					else:
						if mess[0][0][1] == 36:
							self.params["noise_exp"] *= 0.9

						elif mess[0][0][1] == 37:
							self.params["noise_exp"] *= 1/0.9

						if mess[0][0][1] == 33:
							self.params["noise_t"] *= 0.8

						elif mess[0][0][1] == 34:
							self.params["noise_t"] *= 1/0.8

						elif mess[0][0][1] == 39:
							self.params["delay"] *= 1.1

						elif mess[0][0][1] == 38:
							self.params["delay"] *= 1./1.1
						print(self.params)
						engine.set_parameters(self.params)
				elif mess[0][0][0] == 128:
					key_id = mess[0][0][1] - 42 + 18
					engine.key_off(key_id, samplerrate)

			data = engine.get_audio_stream(chunk, samplerrate)
			if record:
				self.audio.append(data)

			max_audio = 10*341193.25
			data *= 0.5/max_audio
			output_bytes = data.tobytes()
			self.stream.write(output_bytes)

			count += 1

# params = {"K": 1., "dynamics": 3.0, "noise_level": 300000., "volume": 600000, "delay": 0.02, "noise_t": 1.0, "Q_mult": 1.0, "t": 0.}
pick_params = {'K': 1.0, 'dynamics': 3.0, 'noise_level': 27144409.1796875, 'volume': 600000,'delay': 0.011289478601075546,'noise_t': 0.2621440000000001,'Q_mult': 1.0,'noise_exp': 1} #pick
fing_params = {'K': 1.0, 'dynamics': 3.0, 'noise_level': 4244409.1796875, 'volume': 600000, 'delay': 0.010, 'noise_t': 1.5500000000001, 'Q_mult': 1, 'noise_exp': 3} #finger



randomize = False

if randomize:
	mult = np.linspace(1., 0., 30)**2*np.random.random(30)

else:
	mult = np.linspace(1., 0., 30)**2

set_multipliers(mult)


a = AudioFile(fing_params)
a.play(record = True, play_midi = True)
engine.reset_timer()








