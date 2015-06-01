/**********************************************************************

  Audacity: A Digital Audio Editor

  ToneGen.cpp

  Steve Jolly
  James Crook (Adapted for 'Chirps')

  This class implements a tone generator effect.

*******************************************************************//**

\class EffectToneGen
\brief An Effect that can generate a sine, square or sawtooth wave.
An extended mode of EffectToneGen supports 'chirps' where the
frequency changes smoothly during the tone.

*//*******************************************************************/

#include "../Audacity.h"

#include <math.h>
#include <float.h>

#include <wx/intl.h>
#include <wx/valgen.h>

#include "../Project.h"
#include "../widgets/NumericTextCtrl.h"
#include "../widgets/valnum.h"

#include "ToneGen.h"

enum kInterpolations
{
   kLinear,
   kLogarithmic,
   kNumInterpolations
};

static const wxString kInterStrings[kNumInterpolations] =
{
   XO("Linear"),
   XO("Logarithmic")
};

enum kWaveforms
{
   kSine,
   kSquare,
   kSawtooth,
   kSquareNoAlias,
   kNumWaveforms
};

static const wxString kWaveStrings[kNumWaveforms] =
{
   XO("Sine"),
   XO("Square"),
   XO("Sawtooth"),
   XO("Square, no alias")
};

// Define keys, defaults, minimums, and maximums for the effect parameters
//
//     Name       Type     Key                  Def      Min      Max                     Scale
Param( StartFreq, double,  XO("StartFreq"),     440.0,   1.0,     DBL_MAX,                1  );
Param( EndFreq,   double,  XO("EndFreq"),       1320.0,  1.0,     DBL_MAX,                1  );
Param( StartAmp,  double,  XO("StartAmp"),      0.8,     0.0,     1.0,                    1  );
Param( EndAmp,    double,  XO("EndAmp"),        0.1,     0.0,     1.0,                    1  );
Param( Frequency, double,  XO("Frequency"),     440.0,   1.0,     DBL_MAX,                1  );
Param( Amplitude, double,  XO("Amplitude"),     0.8,     0.0,     1.0,                    1  );
Param( Waveform,  int,     XO("Waveform"),      0,       0,       kNumWaveforms - 1,      1  );
Param( Interp,    int,     XO("Interpolation"), 0,       0,       kNumInterpolations - 1, 1  );

//
// EffectToneGen
//

BEGIN_EVENT_TABLE(EffectToneGen, wxEvtHandler)
    EVT_TEXT(wxID_ANY, EffectToneGen::OnControlUpdate)
END_EVENT_TABLE();

EffectToneGen::EffectToneGen(bool isChirp)
{
   wxASSERT(kNumWaveforms == WXSIZEOF(kWaveStrings));
   wxASSERT(kNumInterpolations == WXSIZEOF(kInterStrings));

   mChirp = isChirp;

   mWaveform = DEF_Waveform;
   mFrequency[0] = DEF_StartFreq;
   mFrequency[1] = DEF_EndFreq;
   mAmplitude[0] = DEF_StartAmp;
   mAmplitude[1] = DEF_EndAmp;
   mInterpolation = DEF_Interp;

   for (int i = 0; i < kNumWaveforms; i++)
   {
      mWaveforms.Add(wxGetTranslation(kWaveStrings[i]));
   }

   for (int i = 0; i < kNumInterpolations; i++)
   {
      mInterpolations.Add(wxGetTranslation(kInterStrings[i]));
   }
}

EffectToneGen::~EffectToneGen()
{
}

// IdentInterface implementation

wxString EffectToneGen::GetSymbol()
{
   return mChirp
      ? CHIRP_PLUGIN_SYMBOL
      : TONE_PLUGIN_SYMBOL;
}

wxString EffectToneGen::GetDescription()
{
   return mChirp
      ? XO("Generates four different types of tone waveform while allowing starting and ending amplitude and frequency")
      : XO("Generates four different types of tone waveform");
}

// EffectIdentInterface implementation

EffectType EffectToneGen::GetType()
{
   return EffectTypeGenerate;
}

// EffectClientInterface implementation

int EffectToneGen::GetAudioOutCount()
{
   return 1;
}

bool EffectToneGen::ProcessInitialize(sampleCount WXUNUSED(totalLen), ChannelNames WXUNUSED(chanMap))
{
   mPositionInCycles = 0.0;
   mSample = 0;

   fInt = 0.0;

   return true;
}

sampleCount EffectToneGen::ProcessBlock(float **WXUNUSED(inBlock), float **outBlock, sampleCount blockLen)
{
   float *buffer = outBlock[0];
   double throwaway = 0;        //passed to modf but never used
   sampleCount i;
   double f = 0.0;
   double a, b;
   int k;

   double frequencyQuantum;
   double BlendedFrequency;
   double BlendedAmplitude;
   double BlendedLogFrequency = 0.0;

   // calculate delta, and reposition from where we left
   double amplitudeQuantum = (mAmplitude[1] - mAmplitude[0]) / mSampleCnt;
   BlendedAmplitude = mAmplitude[0] + amplitudeQuantum * mSample;

   // precalculations:
   double pre2PI = 2.0 * M_PI;
   double pre4divPI = 4.0 / M_PI;

   // initial setup should calculate deltas
   if (mInterpolation == kLogarithmic)
   {
      // this for log interpolation
      mLogFrequency[0] = log10(mFrequency[0]);
      mLogFrequency[1] = log10(mFrequency[1]);
      // calculate delta, and reposition from where we left
      frequencyQuantum = (mLogFrequency[1] - mLogFrequency[0]) / mSampleCnt;
      BlendedLogFrequency = mLogFrequency[0] + frequencyQuantum * mSample;
      BlendedFrequency = pow(10.0, BlendedLogFrequency);
   }
   else
   {
      // this for regular case, linear interpolation
      frequencyQuantum = (mFrequency[1] - mFrequency[0]) / mSampleCnt;
      BlendedFrequency = mFrequency[0] + frequencyQuantum * mSample;
   }

   // synth loop
   for (i = 0; i < blockLen; i++)
   {
      switch (mWaveform)
      {
      case kSine:
         f = sin(pre2PI * mPositionInCycles / mSampleRate);
         break;
      case kSquare:
         f = (modf(mPositionInCycles / mSampleRate, &throwaway) < 0.5) ? 1.0 : -1.0;
         break;
      case kSawtooth:
         f = (2.0 * modf(mPositionInCycles / mSampleRate + 0.5, &throwaway)) - 1.0;
         break;
      case kSquareNoAlias:
//#define NEWSQUARENOALIAS 1
#ifdef NEWSQUARENOALIAS
         // new method which is more complete at low frequencies, and less efficient at high frequencies.  We need a combination of the two.
         int N = floor((mSampleRate/2)/BlendedFrequency);
         double norm = sqrt((2./N)*BlendedFrequency/mSampleRate);
         double pit = M_PI * mPositionInCycles / mSampleRate;
         double w = cos(N * pit);
         double v = sin((N+1) * pit);
         double u = sin(pit);
         double x, temp;
         if( !(N&1) )   // N even
            if( fabs(u) < 0.000001 )
               x = N+1;
            else
               x = v/u;
         else
            if( fabs(u) < 0.000001 )
               x = (N+1)*cos(pit);
            else
               x = v/u;
         fInt += (w*x - 0.0)*norm;  // there may be a drift in the integration, this needs to be investigated

         double delta = mSampleRate/2.;
         double pit1 = M_PI * (mPositionInCycles + delta) / mSampleRate;
         double w1 = cos(N * pit1);
         double v1 = sin((N+1) * pit1);
         double u1 = sin(pit1);
         if( !(N&1) )   // N even
            if( fabs(u1) < 0.000001 )
               x = N+1;
            else
               x = v1/u1;
         else
            if( fabs(u1) < 0.000001 )
               x = (N+1)*cos(pit1);
            else
               x = v1/u1;
         fInt -= (w1*x - 0.0)*norm;
         f = fInt * 2 - 1.;
#else
         b = (1.0 + cos((pre2PI * BlendedFrequency) / mSampleRate)) / pre4divPI;  //scaling factor for all except fundamental
         //do fundamental (k=1) outside loop
         f = pre4divPI * sin(pre2PI * mPositionInCycles / mSampleRate);
         for (k = 3; (k < 200) && (k * BlendedFrequency < mSampleRate / 2.0); k += 2)
         {
            //Hanning Window in freq domain
            a = 1.0 + cos((pre2PI * k * BlendedFrequency) / mSampleRate);
            if(i==0)
               wxLogDebug(wxT("<point f=\"%f\" d=\"%f\"/>"),k*BlendedFrequency,20.*log10(a));
            //calc harmonic, apply window, scale to amplitude of fundamental
            f += a * sin(pre2PI * mPositionInCycles / mSampleRate * k) / (b * k);
         }
#endif
      }
      // insert value in buffer
      buffer[i] = (float) (BlendedAmplitude * f);
      // update freq,amplitude
      mPositionInCycles += BlendedFrequency;
      BlendedAmplitude += amplitudeQuantum;
      if (mInterpolation == kLogarithmic)
      {
         BlendedLogFrequency += frequencyQuantum;
         BlendedFrequency = pow(10.0, BlendedLogFrequency);
      }
      else
      {
         BlendedFrequency += frequencyQuantum;
      }
   }

   // update external placeholder
   mSample += blockLen;

   return blockLen;
}

bool EffectToneGen::GetAutomationParameters(EffectAutomationParameters & parms)
{
   if (mChirp)
   {
      parms.Write(KEY_StartFreq, mFrequency[0]);
      parms.Write(KEY_EndFreq, mFrequency[1]);
      parms.Write(KEY_StartAmp, mAmplitude[0]);
      parms.Write(KEY_EndAmp, mAmplitude[1]);
   }
   else
   {
      parms.Write(KEY_Frequency, mFrequency[0]);
      parms.Write(KEY_Amplitude, mAmplitude[0]);
   }

   parms.Write(KEY_Waveform, mWaveforms[mWaveform]);
   parms.Write(KEY_Interp, mInterpolations[mInterpolation]);

   return true;
}

bool EffectToneGen::SetAutomationParameters(EffectAutomationParameters & parms)
{
   ReadAndVerifyEnum(Waveform, mWaveforms);
   ReadAndVerifyEnum(Interp, mInterpolations);
   if (mChirp)
   {
      ReadAndVerifyDouble(StartFreq);
      ReadAndVerifyDouble(EndFreq);
      ReadAndVerifyDouble(StartAmp);
      ReadAndVerifyDouble(EndAmp);
      mFrequency[0] = StartFreq;
      mFrequency[1] = EndFreq;
      mAmplitude[0] = StartAmp;
      mAmplitude[1] = EndAmp;
   }
   else
   {
      ReadAndVerifyDouble(Frequency);
      ReadAndVerifyDouble(Amplitude);
      mFrequency[0] = Frequency;
      mFrequency[1] = Frequency;
      mAmplitude[0] = Amplitude;
      mAmplitude[1] = Amplitude;
   }

   mWaveform = Waveform;
   mInterpolation = Interp;

   double freqMax = (GetActiveProject() ? GetActiveProject()->GetRate() : 44100.0) / 2.0;
   mFrequency[1] = TrapDouble(mFrequency[1], MIN_EndFreq, freqMax);

   return true;
}

// Effect implementation

void EffectToneGen::PopulateOrExchange(ShuttleGui & S)
{
   wxTextCtrl *t;

   S.StartMultiColumn(2, wxCENTER);
   {
      wxChoice *c = S.AddChoice(_("Waveform:"), wxT(""), &mWaveforms);
      c->SetValidator(wxGenericValidator(&mWaveform));

      if (mChirp)
      {
         S.AddFixedText(wxT(""));
         S.StartHorizontalLay(wxEXPAND);
         {
            S.StartHorizontalLay(wxLEFT, 50);
            {
               S.AddTitle(_("Start"));
            }
            S.EndHorizontalLay();

            S.StartHorizontalLay(wxLEFT, 50);
            {
               S.AddTitle(_("End"));
            }
            S.EndHorizontalLay();
         }
         S.EndHorizontalLay();

         S.AddPrompt(_("Frequency (Hz):"));
         S.StartHorizontalLay(wxEXPAND);
         {
            S.StartHorizontalLay(wxLEFT, 50);
            {
               FloatingPointValidator<double> vldStartFreq(6, &mFrequency[0], NUM_VAL_NO_TRAILING_ZEROES);
               vldStartFreq.SetRange(MIN_StartFreq, GetActiveProject()->GetRate() / 2.0);
               t = S.AddTextBox(wxT(""), wxT(""), 12);
               t->SetName(_("Frequency Hertz Start"));
               t->SetValidator(vldStartFreq);
            }
            S.EndHorizontalLay();

            S.StartHorizontalLay(wxLEFT, 50);
            {
               FloatingPointValidator<double> vldEndFreq(6, &mFrequency[1], NUM_VAL_NO_TRAILING_ZEROES);
               vldEndFreq.SetRange(MIN_EndFreq, GetActiveProject()->GetRate() / 2.0);
               t = S.AddTextBox(wxT(""), wxT(""), 12);
               t->SetName(_("Frequency Hertz End"));
               t->SetValidator(vldEndFreq);
            }
            S.EndHorizontalLay();
         }
         S.EndHorizontalLay();

         S.AddPrompt(_("Amplitude (Hz):"));
         S.StartHorizontalLay(wxEXPAND);
         {
            S.StartHorizontalLay(wxLEFT, 50);
            {
               FloatingPointValidator<double> vldStartAmp(6, &mAmplitude[0], NUM_VAL_NO_TRAILING_ZEROES);
               vldStartAmp.SetRange(MIN_StartAmp, MAX_StartAmp);
               t = S.AddTextBox(wxT(""), wxT(""), 12);
               t->SetName(_("Amplitude Start"));
               t->SetValidator(vldStartAmp);
            }
            S.EndHorizontalLay();

            S.StartHorizontalLay(wxLEFT, 50);
            {
               FloatingPointValidator<double> vldEndAmp(6, &mAmplitude[1], NUM_VAL_NO_TRAILING_ZEROES);
               vldEndAmp.SetRange(MIN_EndAmp, MAX_EndAmp);
               t = S.AddTextBox(wxT(""), wxT(""), 12);
               t->SetName(_("Amplitude End"));
               t->SetValidator(vldEndAmp);
            }
            S.EndHorizontalLay();
         }
         S.EndHorizontalLay();

         c = S.AddChoice(_("Interpolation:"), wxT(""), &mInterpolations);
         c->SetValidator(wxGenericValidator(&mInterpolation));
      }
      else
      {
         FloatingPointValidator<double> vldFrequency(6, &mFrequency[0], NUM_VAL_NO_TRAILING_ZEROES);
         vldFrequency.SetRange(MIN_Frequency, GetActiveProject()->GetRate() / 2.0);
         t = S.AddTextBox(_("Frequency (Hz):"), wxT(""), 12);
         t->SetValidator(vldFrequency);

         FloatingPointValidator<double> vldAmplitude(6, &mAmplitude[0], NUM_VAL_NO_TRAILING_ZEROES);
         vldAmplitude.SetRange(MIN_Amplitude, MAX_Amplitude);
         t = S.AddTextBox(_("Amplitude (0-1):"), wxT(""), 12);
         t->SetValidator(vldAmplitude);
      }

      bool isSelection;
      double duration = GetDuration(&isSelection);

      S.AddPrompt(_("Duration:"));
      mToneDurationT = new
         NumericTextCtrl(NumericConverter::TIME,
                        S.GetParent(),
                        wxID_ANY,
                        isSelection ? _("hh:mm:ss + samples") : _("hh:mm:ss + milliseconds"),
                        duration,
                        mProjectRate,
                        wxDefaultPosition,
                        wxDefaultSize,
                        true);
      mToneDurationT->SetName(_("Duration"));
      mToneDurationT->EnableMenu();
      S.AddWindow(mToneDurationT, wxALIGN_LEFT | wxALL);
   }
   S.EndMultiColumn();

   return;
}

bool EffectToneGen::TransferDataToWindow()
{
   if (!mUIParent->TransferDataToWindow())
   {
      return false;
   }

   mToneDurationT->SetValue(GetDuration());

   return true;
}

bool EffectToneGen::TransferDataFromWindow()
{
   if (!mUIParent->Validate() || !mUIParent->TransferDataFromWindow())
   {
      return false;
   }

   if (!mChirp)
   {
      mFrequency[1] = mFrequency[0];
      mAmplitude[1] = mAmplitude[0];
   }

   SetDuration(mToneDurationT->GetValue());

   return true;
}

// EffectToneGen implementation

void EffectToneGen::OnControlUpdate(wxCommandEvent & WXUNUSED(evt))
{
   if (!EnableApply(mUIParent->TransferDataFromWindow()))
   {
      return;
   }
}
