using SeisProcessing
using Test


d = SeisLinearEvents();
snr_i =1.2
dn = SeisAddNoise(d,snr_i);

snr = MeasureSNR(d,dn)

# test that snr = snr_i
quality_factor = abs(snr-snr_i)
@test quality_factor < 1