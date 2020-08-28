#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include "waveforms.c"

#define ARRAY_COUNT(a) (sizeof(a) / sizeof(a[0]))

void OutputSignal(const char *filename, double *signal, int count) {
    FILE *OutputFile = fopen(filename, "w");
    assert(OutputFile != NULL);
    for (int i = 0; i < count; ++i) {
        fprintf(OutputFile, "\n%f", signal[i]);
    }
    fclose(OutputFile);
}

double SignalMean(double *signal, int count){
    double mean = 0.0;
    
    for (int i = 0; i < ARRAY_COUNT(InputSignal_f32_1kHz_15kHz); ++i) {
        mean += signal[i];
    }
    mean = mean / count;
    return mean;
}

double SignalVariance(double mean, double *signal, int count) {
    
    double variance = 0.0;
    for (int i = 0; i < count; ++i) {
        double val = signal[i] - mean;
        variance += val * val;
    }
    
    variance /= (count - 1);
    return variance;
}

double StandardDeviation(double variance) {
    return sqrt(variance);
}

void RunningSum(double *input, int input_size, double *output) {
    output [0] = input [0];
    for (int i = 1; i < input_size; ++i) {
        output[i] = input[i] + output[i-1];
    }
}

void Convolution(double *Input, int InputSize, double *Impulse, int ImpulseSize, double *Output) {
    for (int i = 0; i < InputSize; ++i) {
        for (int j = 0; j < ImpulseSize; ++j) {
            Output[i + j] += Input[i] * Impulse[j];
        }
    }
}

void DFT(double *input, int input_size, double *output_rel, double *output_im) {
    int i_max = input_size, k_max = input_size / 2;
    
    for (int k = 0; k < k_max; ++k) {
        output_rel[k] = 0;
        output_im[k] = 0;
    }
    
    double one_over_n = 1.0 / ((double)input_size);
    
    for (int k = 0; k < k_max; ++k) {
        for (int i = 0; i < i_max; ++i) {
            double tpki = 2.0 * M_PI * k * i;
            output_rel[k] = output_rel[k] + input[i] * cos(tpki * one_over_n);
            output_im[k]  = output_im[k] - input[i] * sin(tpki * one_over_n);
        }
    }
}

void IDFT(double *input_rel, double *input_im, int input_size, double *output) {
    int n = input_size * 2, n2 = input_size;
    double one_over_n2 = 1.0 / ((double)n2);
    double one_over_n = 1.0/ ((double)n);
    
    for (int i = 0; i < n; ++i) output[i] = 0;
    
    double *input_brel = malloc(n2 * sizeof(double));
    double *input_bim = malloc(n2 * sizeof(double));
    
    input_brel[0] = input_rel[0] / n;
    input_bim[0] = -input_im[0] / n;
    
    for (int k = 1; k < n2; ++k) {
        input_brel[k] = input_rel[k] * one_over_n2;
        input_bim[k] = -input_im[k] * one_over_n2;
    }
    
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n2; ++k) {
            double tpki = 2.0 * M_PI * k * i;
            output[i] += input_brel[k] * cos(tpki * one_over_n) + input_bim[k] * sin(tpki * one_over_n);
        }
    }
    
    free(input_brel);
    free(input_bim);
}

void DFTMagnitude(const double *Rel, const double *Im, int input_size, double *output) {
    
    for (int k = 0; k < input_size; ++k) output[k] = 0;
    
    for (int k = 0; k < input_size; ++k) {
        double r = Rel[k];
        double i = Im[k];
        output[k] = sqrt(r*r + i*i);
    }
}

void RectToPolar(double *Rel, double *Im, int InputSize, double *Mag, double *Phase) {
    for (int i = 0; i < InputSize; ++i) {
        Mag[i] = sqrt(Rel[i] * Rel[i] + Im[i] * Im[i]);
        
        if (fabs(Rel[i]) < 0.0000001) Rel[i] = 0.0000000000000000000000000001;
        Phase[i] = atan(Im[i] / Rel[i]);
        
        if (Rel[i] < 0.0 && Im[i] < 0.0) Phase[i] -= M_PI;
        if (Rel[i] < 0.0 && Im[i] >= 0.0) Phase[i] += M_PI;
    }
}

void PolarToRect(const double *Mag, const double *Phase, int InputSize, double *Rel, double *Im) {
    for (int i = 0; i < InputSize; ++i) {
        Rel[i] = Mag[i] * cos(Phase[i]);
        Im[i] = Mag[i] * sin(Phase[i]);
    }
}

void ComplexDFT(const double *input_td_rex, const double *input_td_imx, int input_size, double *output_fd_rex, double *output_fd_imx) {
    
    double SR, SI;
    double inv_len = (1.0 / (double)input_size);
    for (int k = 0; k < input_size - 1; ++k) {
        for (int i = 0; i < input_size - 1; ++i) {
            SR = cos(2.0 * M_PI * (double)k * (double)i * inv_len);
            SI = -sin(2.0 * M_PI * (double)k * (double)i * inv_len);
            
            output_fd_rex[k] += input_td_rex[i] * SR - input_td_imx[i] * SI;
            output_fd_imx[k] += input_td_imx[i] * SI - input_td_rex[i] * SR;
        }
    }
}

//Blackman Window
//Hammin Window

static inline double HammingWindow(double ratio) {
    return 0.54 - 0.46 * cos(2.0 * M_PI * ratio);
}

void LowpassWindowedSincFilterKernel(double Cutoff,
                                     int KernelSize,
                                     double *OutKernel) {

    assert(KernelSize & 1); // Should be odd
    int Window = KernelSize - 1;
    int Center = Window / 2;
    double x = 2.0 * M_PI * Cutoff;
    for (int i = 0; i < KernelSize; ++i) {
        int Distance = i - Center;
        if (Distance == 0 ) OutKernel[i] = x;
        else {
            OutKernel[i] = sin(x * ((double)Distance)) / ((double)Distance);
            OutKernel[i] *= HammingWindow((double)i / (double)Window); // HammingWindow
        }
    }
}

void SpectralInversion(double *Kernel, int KernelSize) {
    for (int i = 0; i < KernelSize; ++i) {
        Kernel[i] = -(Kernel[i]);
    }
    Kernel[KernelSize / 2] += 1.0;
}

void SpectralReversal(double *Kernel, int KernelSize) { 
    double v = 1.0;
    for (int i = 0; i < KernelSize; ++i) {
        Kernel[i] *= v;
        v *= -1.0;
    }
}

void BandpassWindowedSincFilterKernel( double *LowCutoffBuffer, double *HighCutoffBuffer, double *Kernel, int KernelSize, double LowCutoff, double HighCutoff) {
    
    LowpassWindowedSincFilterKernel(LowCutoff, KernelSize, LowCutoffBuffer);
    LowpassWindowedSincFilterKernel(HighCutoff, KernelSize, HighCutoffBuffer);
    SpectralInversion(HighCutoffBuffer, KernelSize);
    
    for (int i = 0; i < KernelSize; ++i) {
        Kernel[i] = LowCutoffBuffer[i] + HighCutoffBuffer[i];
    }
    
    SpectralInversion(Kernel, KernelSize);
}

int main(int argc, char **argv) {
    (void)argc; (void)argv;
    
#if 0
    // Simple statistics
    double mean = SignalMean(&InputSignal_f32_1kHz_15kHz[0], ARRAY_COUNT(InputSignal_f32_1kHz_15kHz));
    printf("Mean: %f\n", mean);
    
    double variance = SignalVariance(mean, &InputSignal_f32_1kHz_15kHz[0], ARRAY_COUNT(InputSignal_f32_1kHz_15kHz));
    printf("Variance: %f\n", variance);
    
    double standard_deviation = StandardDeviation(variance);
    printf("Standard Deviation: %f\n", standard_deviation);
#endif
    
#if 0   
    
#define SignalLength (ARRAY_COUNT(InputSignal_f32_1kHz_15kHz))
#define ImpulseLength (ARRAY_COUNT(Impulse_response))
    OutputSignal("output/input_signal.dat", &InputSignal_f32_1kHz_15kHz[0],SignalLength);
    OutputSignal("output/impulse_response.dat", &Impulse_response[0],ImpulseLength);
    
    double 
        ConvolvedSignal[SignalLength + ImpulseLength] = {0};
    
    Convolution(&InputSignal_f32_1kHz_15kHz[0], SignalLength, &Impulse_response[0], ImpulseLength, &ConvolvedSignal[0]);
    OutputSignal("output/output_signal.dat", &ConvolvedSignal[0],ARRAY_COUNT(ConvolvedSignal));
    
#endif
    
#if 0    
    double running_sum[ARRAY_COUNT(InputSignal_f32_1kHz_15kHz)] = {0};
    RunningSum(InputSignal_f32_1kHz_15kHz, ARRAY_COUNT(InputSignal_f32_1kHz_15kHz), &running_sum[0]);
    OutputSignal("output/output_runningsum.dat", &running_sum[0],ARRAY_COUNT(running_sum));
#endif
    
#if 0    
    double RelX[ARRAY_COUNT(InputSignal_f32_1kHz_15kHz) / 2] = {0}; 
    double ImX[ARRAY_COUNT(RelX)] = {0};
    double Mag[ARRAY_COUNT(ImX)] = {0};
    DFT(&InputSignal_f32_1kHz_15kHz[0], ARRAY_COUNT(InputSignal_f32_1kHz_15kHz), &RelX[0], &ImX[0]);
    DFTMagnitude(RelX, ImX, ARRAY_COUNT(RelX), Mag);
    OutputSignal("output/input_signal.dat", &InputSignal_f32_1kHz_15kHz[0], ARRAY_COUNT(InputSignal_f32_1kHz_15kHz));
    OutputSignal("output/output_rel.dat", &RelX[0], ARRAY_COUNT(RelX));
    OutputSignal("output/output_im.dat", &ImX[0], ARRAY_COUNT(ImX));
    OutputSignal("output/output_mag.dat", &Mag[0], ARRAY_COUNT(Mag));
    
    double OutputIDFT[ARRAY_COUNT(InputSignal_f32_1kHz_15kHz)] = {0};
    IDFT(&RelX[0], &ImX[0], ARRAY_COUNT(RelX), &OutputIDFT[0]);
    
    OutputSignal("output/output_idft.dat", &OutputIDFT[0], ARRAY_COUNT(OutputIDFT));
#endif
    
#if 0
    double RelX[ARRAY_COUNT(_640_points_ecg_) / 2] = {0}; 
    double ImX[ARRAY_COUNT(_640_points_ecg_) / 2] = {0};
    
    DFT(&_640_points_ecg_[0], ARRAY_COUNT(_640_points_ecg_), &RelX[0], &ImX[0]);
    
    OutputSignal("output/input_ecg.dat", &_640_points_ecg_[0], ARRAY_COUNT(_640_points_ecg_));
    OutputSignal("output/output_rel.dat", &RelX[0], ARRAY_COUNT(RelX));
    OutputSignal("output/output_im.dat", &ImX[0], ARRAY_COUNT(ImX));
    
    double OutputIDFT[ARRAY_COUNT(_640_points_ecg_)] = {0};
    IDFT(&RelX[0], &ImX[0], ARRAY_COUNT(RelX), &OutputIDFT[0]);
    
    OutputSignal("output/output_idft.dat", &OutputIDFT[0], ARRAY_COUNT(OutputIDFT));
#endif
    
#if 0    
#define SIG_LEN (ARRAY_COUNT(_320_pts_ecg_REX))
    
    double Mag[SIG_LEN] = {0};
    double Phase[SIG_LEN] = {0};
    
    RectToPolar(&_320_pts_ecg_REX[0],&_320_pts_ecg_IMX[0], SIG_LEN, &Mag[0], &Phase[0]);
    
    OutputSignal("output/output_320ecg_rex.dat", &_320_pts_ecg_REX[0], SIG_LEN);
    OutputSignal("output/output_320ecg_imx.dat", &_320_pts_ecg_IMX[0], SIG_LEN);
    OutputSignal("output/output_320ecg_Mag.dat", &Mag[0], SIG_LEN);
    OutputSignal("output/output_320ecg_Phase.dat", &Phase[0], SIG_LEN);
#endif
    
#if 0    
#define SIG_LEN (ARRAY_COUNT(_501pts_20Hz_sig_rex))
    
    double OutREX[SIG_LEN] = {0};
    double OutIMX[SIG_LEN] = {0};
    
    ComplexDFT(&_501pts_20Hz_sig_rex[0], &_501pts_20Hz_sig_imx[0], SIG_LEN, &OutREX[0], &OutIMX[0]);
    
    OutputSignal("output/input_501pts_rex.dat", &_501pts_20Hz_sig_rex[0], SIG_LEN);
    OutputSignal("output/input_501pts_imx.dat", &_501pts_20Hz_sig_imx[0], SIG_LEN);
    OutputSignal("output/output_501pts_rex.dat", &OutREX[0], SIG_LEN);
    OutputSignal("output/output_501pts_imx.dat", &OutIMX[0], SIG_LEN);
#endif
    
#define SignalLength (ARRAY_COUNT(InputSignal_f32_1kHz_15kHz))
#define KernelSize 29
#define OutputSize (SignalLength)
#define Input InputSignal_f32_1kHz_15kHz

    double Cutoff = 0.1;
    
    double LowpassKernel[KernelSize] = {0};
    double OutputLowpass[OutputSize] = {0};
    
    LowpassWindowedSincFilterKernel(Cutoff, KernelSize, LowpassKernel);
    SpectralInversion(LowpassKernel, KernelSize);
    for (int i = 0; i < SignalLength; ++i) {
        OutputLowpass[i] = 0;
        for (int j = 0; j < KernelSize; ++j) {
            OutputLowpass[i] += Input[i - j] * LowpassKernel[j];
        }
    }
    OutputSignal("output/output_filter_kernel_lowpass.dat", &LowpassKernel[0], KernelSize);
    OutputSignal("output/output_filter_signal_lowpass.dat", &OutputLowpass[0], SignalLength);

    double FreqResp_REX[KernelSize] = {0};
    double FreqResp_IMX[KernelSize] = {0};
    DFT(LowpassKernel, KernelSize, FreqResp_REX, FreqResp_IMX);

    double Mag[KernelSize] = {0};
    double Phase[KernelSize] = {0};
    RectToPolar(FreqResp_REX, FreqResp_IMX, KernelSize, Mag, Phase);

    OutputSignal("output/output_filter_kernel_mag.dat", &Mag[0], KernelSize);
    OutputSignal("output/output_filter_kernel_phase.dat", &Phase[0], KernelSize);

#if 0
    double BandpassKernel[KernelSize] = {0};
    double OutputBandpass[OutputSize] = {0};
    double Buffer1[KernelSize] = {0};
    double Buffer2[KernelSize] = {0};

    double LowCutoff = 0.1;
    double HighCutoff = 0.2;

    BandpassWindowedSincFilterKernel(Buffer1, Buffer2, BandpassKernel, KernelSize, LowCutoff, HighCutoff);
    for (int i = 0; i < SignalLength; ++i) {
        OutputBandpass[i] = 0;
        for (int j = 0; j < KernelSize; ++j) {
            OutputBandpass[i] += ((i - j) < 0 ? 0.0 : Input[i - j]) * BandpassKernel[j];
        }
    }
    
    OutputSignal("output/output_filter_kernel_bandpass.dat", &BandpassKernel[0], KernelSize);
    OutputSignal("output/output_filter_signal_bandpass.dat", &OutputBandpass[0], SignalLength);
    
    double FreqResp_REX[KernelSize] = {0};
    double FreqResp_IMX[KernelSize] = {0};
    
    DFT(LowpassKernel, KernelSize, FreqResp_REX, FreqResp_IMX);
    OutputSignal("output/output_filter_freqresp_rex.dat", FreqResp_REX, KernelSize);
    OutputSignal("output/output_filter_freqresp_imx.dat", FreqResp_IMX, KernelSize);
#endif
    
    return 0;
}
