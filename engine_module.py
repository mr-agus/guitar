#initialize_harmonic_info
import os
import numpy as np
import time

import numpy.ctypeslib as npct
from ctypes import c_int, c_float

array_1d_float = npct.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
array_1d_int   = npct.ndpointer(dtype=np.int32,   ndim=1, flags='C_CONTIGUOUS')



if not os.path.isdir("lib"):
	os.system("mkdir lib")

os.system("gcc -c -Ofast -fPIC source/engine.c -o lib/engine.o -lm -fopenmp -std=c99")
os.system("gcc -shared lib/engine.o -o lib/engine.so -lm -fopenmp -std=c99")


libcd = npct.load_library("lib/engine.so", ".")

deg2rad = np.pi/180

libcd.load_harmonic_multipliers.restype = None
libcd.load_harmonic_multipliers.argtypes = [array_1d_float, c_int]

libcd.initialize_harmonic_info.restype = None
libcd.initialize_harmonic_info.argtypes = [c_int, c_int, array_1d_float, array_1d_float, array_1d_float, array_1d_float, array_1d_float, array_1d_float, array_1d_float, array_1d_float, array_1d_float, array_1d_float, array_1d_int]

libcd.fill_chunk.restype = None
libcd.fill_chunk.argtypes = [array_1d_float, c_int, c_int]

libcd.play_note.restype = None
libcd.play_note.argtypes = [c_int, c_float, c_int]

libcd.key_off.restype = None
libcd.key_off.argtypes = [c_int, c_int]




libcd.initialize_oscillators.restype = None
libcd.initialize_oscillators.argtypes = []

libcd.reset_timer.restype  = None
libcd.reset_timer.argtypes = []

libcd.set_parameters.restype = None
libcd.set_parameters.argtypes = [c_float, c_float, c_float, c_float, c_float, c_float, c_float, c_float]
#initialize_harmonic_info(int min_key_value, unsigned int key_amount, float *freqs, float *t_decays, float *amplitudes_left, float *amplitudes_right, float *b1s, float *b2s, float *ks, float *f1s, float *values_l, float *values_r)

#set_parameters(float dynamics, float noise_level, float volume, float K, float delay, float noise_t, float Q_mult)

#play_note(unsigned int key_id, float amplitude)

def reset_timer():
	libcd.reset_timer()

def set_parameters(params):
	dynamics, noise_level, volume, K, delay, noise_t, Q_mult, noise_exp = params["dynamics"], params["noise_level"], params["volume"], params["K"], params["delay"], params["noise_t"], params["Q_mult"], params["noise_exp"]
	libcd.set_parameters(np.float32(dynamics), np.float32(noise_level), np.float32(volume), np.float32(K), np.float32(delay), np.float32(noise_t), np.float32(Q_mult), np.float32(noise_exp))

def set_multipliers(multipliers):
	multipliers = np.ascontiguousarray(multipliers, dtype = np.float32)
	libcd.load_harmonic_multipliers(multipliers, len(multipliers))

def play_note(key_id, amplitude, samplerrate):
	libcd.play_note(np.int32(key_id), np.float32(amplitude), np.int32(samplerrate))

def key_off(key_id, samplerrate):
	libcd.key_off(np.int32(key_id), np.int32(samplerrate))

def initialize_engine():
	libcd.initialize_oscillators()

def load_harmonics(min_key_value, key_amount, freqs, t_decays, amplitudes_left, amplitudes_right, b1s, b2s, ks, f1s, values_l, values_r, nharms):
	min_key_value    = np.int32(min_key_value)
	key_amount       = np.int32(key_amount)

	freqs            = np.ascontiguousarray(freqs,            dtype = np.float32)
	t_decays         = np.ascontiguousarray(t_decays,         dtype = np.float32)
	amplitudes_left  = np.ascontiguousarray(amplitudes_left,  dtype = np.float32)
	amplitudes_right = np.ascontiguousarray(amplitudes_right, dtype = np.float32)
	b1s              = np.ascontiguousarray(b1s,              dtype = np.float32)
	b2s              = np.ascontiguousarray(b2s,              dtype = np.float32)
	ks               = np.ascontiguousarray(ks,               dtype = np.float32)
	f1s              = np.ascontiguousarray(f1s,              dtype = np.float32)
	values_l         = np.ascontiguousarray(values_l,         dtype = np.float32)
	values_r         = np.ascontiguousarray(values_r,         dtype = np.float32)
	nharms           = np.ascontiguousarray(nharms,           dtype = np.int32)

	libcd.initialize_harmonic_info(min_key_value, key_amount, freqs, t_decays, amplitudes_left, amplitudes_right, b1s, b2s, ks, f1s, values_l, values_r, nharms)

def get_audio_stream(chunk_size, samplerrate):
	results = np.zeros(chunk_size*2, dtype = np.float32)

	libcd.fill_chunk(results, np.int32(chunk_size), np.int32(samplerrate))

	return np.reshape(results, [chunk_size, 2])

