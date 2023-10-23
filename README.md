# The Music Demixing Machine --  NTF-panning

Multichannel music source separation based on NMF (Non-Negative Matrix Factorization) with panning information. The method is intended for classical music recorded with close microphones. This is the MATLAB implementation of the system described in this [paper](https://link.springer.com/article/10.1007/s11227-023-05192-5).

## Model initialization

Run the function `main_train(INITFOLDER)` to initialize the parameters of the model (panning matrix *M_js* and bases *S_pfj*). The input directory `INITFOLDER` must contain:

- WAV files (at 44.1 kHz) with a multichannel recording (one WAV per microphone) contaning instrument solo segments.
- A `marker.txt` annotation file indicating the instrument codes and time marks where each instrument starts playing solo.

### Marker file

Here it is an example of `marker.txt`:

```
00:00:00:000	TR_I
00:00:10:000	HR
00:00:20:000	TB
00:00:30:000	TR_II
```

The time marks are written in the format *hh:mm:ss:ms*. In this example, the trumpet *TR_I* plays solo from 0 to 10 seconds, and then the french horn *HR* from 10 to 20, etc. These segments are used to learn the panning coefficients for each instrument.

You can also specify segments with no solos, like this:

```
00:00:00:000	TR_I
00:00:05:000	SILENCE
00:00:10:000	HR
00:00:20:000	TB
00:00:30:000	TR_II
```

Here, the `main_train` function will just ignore the interval from 5 to 10 seconds .

See the folder `INSTRUMENT_BASES ` to see all the available two-letter instrument codes. Observe than you can add a suffix to the instrument code to identify different sections/performers (*TR_I*, *TR_II*, *TR_A*, *TR_B*, ...).


## Separation

To run the separation, simply call the function `main_separate(INFOLDER, INITFOLDER, OUTFOLDER)`.

- `INFOLDER`: input folder with the multichannel mixture to separate  (WAV files, 44.1 kHz).
- `INITFOLDER`: folder with model parameters (the same used in the previous section).
- `OUTFOLDER`: output folder.


## Execution Example

An out-of-the box example recording is available in the `URMP` folder for testing the method. First, initialize the model:

```bash
initfolder = './URMP/31_Slavonic_tpt_tpt_hn_tbn/train';
main_train('./URMP/31_Slavonic_tpt_tpt_hn_tbn/train');
```

Then, separate the instruments:

```bash
infolder = './URMP/31_Slavonic_tpt_tpt_hn_tbn/mixtures';
outfolder = './OUTPUT';

main_separate(infolder, initfolder, outfolder);
```

