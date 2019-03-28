# He* Spectroscopy

## 2019-03-18
### Data acquisition
* Connected MOT photodiode to Analog in, now records via USBEDev1/ai2, replacing the compressed SFP scan
	* Impedance mismatch between PD/CRO and DAQ - when connected to DAQ, PD trace is v noisy
		* Tried permuting connections & adding terminator, no major improvement
		* Signal acquired by analog in seems relatively clean; has probably been LPF'd
		* The noise manifests spikes at T~45us, with smaller ones at T~10us, but not really visible in analog in. So, not really important now but should check out in future


### Hardware
* Tuned laser to  363.653317MHz (red) - tried to couple blue back in but can't seem to get light through fibre, even with the alignment tool?! Suspect possibly not spec'd for red
	* Will try later with 432nm 'spotlight'
* Removed polarization optics from insertion window

### Pulse sequence

Currently, MOT flashes to resonance at 0.4686 seconds in the 0.6s transfer.
We would prefer a longer hold to allow for greater interrogation time. 


Old 0.6s transfersettings201813Nov201731.xml
Note in old sequence, quad coils switch off 100us before resonant pulse

Extend MOT hold out to 900ms
settings201918Mar141550

c:\remote\settings201918Mar143009.xml

100ms extended 0.6s settings201918Mar144419
1sec load settings201918Mar145542

NB in normal BEC settings 6289 ch 0, 1 ,2 , 7, 8 switch at end of MOT load

Probe beam off c:\remote\settings201918Mar155231.xml

c:\remote\settings201919Mar161045.xml

## 2019-03-19

Calculating scan parameters;

Centre wavelength 412.0815 in air
363.644749 THz in vacuum
Prediction given to +-100fm


Had a quick go, scanning +-250MHz (red) in 2MHz steps, nothing obvious in scope. Will set up a wider scan to run overnight:

NB restructured the interface code, now configured by changing # meas and # cal shots, then computes update interval appropriately

Laser power setpoint struggling; output power from AO is 220mW, with 90mW at fibre input and 24mW at fibre exit. Pretty poor! Can worry about efficiency later (maybe need to polish fibre...), but definitely need to fix up set point. BRB, going for dinner.

Todo
* Fix set point DONE
	* Switched PID and RF Attenuator to new PSU. Tried running separately but PID appears bunk on exiting PSU (now on test bench). Should test later; leftmost unit was failing to produce voltage when in run mode, even when all other devices disconnected. Currently running below PID spec voltage to stay within spec of the RF attenuator, but loop appears to be working ok. Outstanding issue: Terrible coupling efficiency
* coupling efficiency FIXED
	* Performed some open-heart surgery on the cavity while measuring power at fibre input. Looks like beam was clipping on the case lid on the way out! Shuffling it a couple mm provided a beautiful beam spot and a nice high power throughput. Some alignment later and we have fibre coupling of 7.5/18 = 42% at low power, up to 80/220 (36%) at high power

Setting up for overnight scan:
* Adjusted power monitor HWP to give steady voltage reading 750-800mW
* Set interrogation setpt 0.5V for interrogation sequence
* Find replacement shutter in optics drawers - burned spot on the face of the shutter, but it works! Now shuts off push & ZS beams at 0.45sec in sequence. For good measure, switch off MOT1 & collimator beams at 0.45s. i.e. ALL loading optics now dark after 0.45 sec, with resonant flash at 0.97sec => 0.52s interrogation time (although probe on during entire sequence)
* average fluro height 90mV
* Probe flipper mirror now set to negative polarity in calibration sequence so it stays inserted during all calibration runs

Overnight scan parameters:
Scanning +- 1GHz, 3MHz steps (both in red), 1 calibration & 4 test shots per wavelength 
(No chance of getting both 412nm transitions out in one scan, they're 30GHz apart..)

=> 5 sec per wavelength
=> 53 minutes per scan. 
10 hours runtime (until 0930) => 11 runs
Also, that's about 35k files to generate!!! ~ 6GB of data

Fiddlesticks: What about the polarization of the light?
We're trying to drive the 3P2 to 3S1 transition, so delta_L = -1, meaning we need left(??) circularly polarized light in the atomic frame! 
However, the magnetic field will be pointing all over the place, and quite weak, so  we might just have to make do with linear light for now and hope that some fraction of the atoms see that as the required polarization; assuming we're even hitting the damn thing.

## 2019-03-20
### Good morning, Kieran!
Hopefully the run was stable overnight. I took a huge scan, as you can see. With luck we'll find something in there. The first thing one could do would be to beef up the analysis script. 
	* Import from a directory with tonight's data
	* Import the WM log and ensure the probe set points can be associated with the analog files
	* Create two structures: one for the calibration data and one for the test data
	* For the calibration data, ensure the time traces are easily aligned and find the average of all calibration shots
	* For the test data, collect all the shots for a given setpoint, align their traces & take an average
	* As a first shot at a metric, compute a few simple things like:
		* Peak height
		* Decay rate of exponential fit to fluro
		* Pointwise mean square error (L^2 norm of the difference) between a given wavelength and the calibration data
		* Plot all of these versus setpoint
		* Crack a cold one/be sad, conditioned on outcome
	* If you get through all this, have Bryce walk you through the wavemeter calibration
	* Also have a think about whether polarization matters here, and what to do about it...
	* Further work that will need to be done for 'real' runs:
		* check/rejig the analog in to ensure we get the probe regulation PD and scanning fabry-perot data. 
		* We'll eventually need a comparison between the WM log and the desired setpoint.
	* I should be in by lunchtime if not before (signing out at midnight)
	* Good hunting!

### JR 
* Beginnings of analysis written; pared back analog import from tuneout and refactored a little bit
* Notice that importing and caching a huge number of files takes aaaaages - thousands of save calls!
	* try adding an option to cache single imports, testing on a smaller directory, importing 1192 files
		* cache_single = true : 
			* 1193 calls to function_cache, 82s total, 53s self time, 83s total execution
			* 522 calls ai_log_single, 8.7s total time, 6 sec sec self time
		* cache_single = false : 
			* 1 call to function_cache, 24s total time, 2.7s self time, 25sec total execution
			* 20 sec ai_log_single, 14.5 sec self time
	* Wait, why the discrepancy in calls to ai_log_single?
		* OH. cache_opts aren't actually passed from ai_log_import to ai_log_import_core; they're configured differently
			* Removing all internal config from functions, this is cancer. This problem took me like two hours
		* Try again. 200 files
			* cache_single_import enabled
				* 11.04 sec in function_cache
			* cache_single_import disabled
				* 3.2 sec total
				* 0.527s function_cache self time
			Verdict: CRUSHED
				* Granted, once the import is called, it shouldn't need to be done again; but this is still way faster
				* Processing will now be done by second-stage function with fixed import, more agile dev speedy speedy good!
			* Testing with the big run: Took tens of thousands of seconds last time, let's give it a go
				* ~ 5min, factor of 6 speed up

* Moving on; Good news is the WM log function is unchanged so that just works. Moved configs out of function and into config file, otherwise unchanged and good to go
* ditto LV log, just ripped it out of tuneout main and made it a function import_lv_log, might be worth generalizing a bit

### Extending MOT hold
avg of 4 shots after running for 10
Hold time relative to +520ms
T (ms)	Fluro peak (mv)
0		111
100		105
300		104
1000	120?!
ugh problems


## Hello there!

https://tinyurl.com/y5wt6nvj

So, it turns out that the sequence I set up the other night was not what it seemed.
Turns out that it was 450ms of darkness, then loading a MOT and fluoro as normal. Oops.
It's important that we fix this; with a closing window for measurements, runtime is precious!
I discovered this by extending the MOT "hold" and noting that the fluro did not decrease, as you see above.
The diagnosis is:
	* The shutter I replaced ALSO switches the horizontal beams on the second MOT, so sadly we can't use it to turn off MOT loading
	* The easiest solution, I think, is just to switch off the MOT1 beams after the load. I suggest free-running the sequence and walking around with an IR card; look at the beam entry ports, look at the exit apertures on the optics table. Get a feel for the timing of everything, then compare your observations to the pulse sequences in LabView. 
	* The goal is to ensure that the sequence is composed of
		* MOT load: All beams on
		* MOT hold: MOT1 beams switch off
		* Fluro: MOT2 beams flash close to resonance, probe beam switching off somehow to avoid reflections 
Will probably involve a lot of trial and error, but that's how things get done around here ;)

So, task for today: Fix the sequence! Some tips:
	* Otherwise and afterwards: Determine the centre frequency and scan range for a finer search around the expected wavelength. Bonus points if you do this for the other target transitions too. 
	* When both these things are done, get the scan going! Then keep pushing the code forward

## (2019-03-21) Moves towards seeing something


Determine MOT lifetime to be approx 250ms, so use that as test sequence

160 ms hold
	Interagation sequence: c:\remote\settings201921Mar110433.xml
	Calibration sequence: c:\remote\settings201921Mar111103.xml

260 ms hold
	Interagation sequence: c:\remote\settings201921Mar113737.xml
	Calibration sequence: c:\remote\settings201921Mar113847.xml


### Calibrating wavemeter
consider using two transitions!

Stop WM feedback
note setpt
tisaf setpt to 822.2 
run wm feedback without blue
run measure_2photon_transitions
turn on pmt
run measure ~ 3 times
stash in folder
stop WM feedback
turn on modulation fn 
select ch1 on settings signal
stop wm monitor
Operation>Calibration>etc0
kill matlab 
run ws8 fb then meas_2p 
Kill meas_2p, reboot ws8 feedback & restore noted setpt
power off pmt & mod

### Night works
Having trouble keeping laser locked :(
Yeah, spent a few hours battling the laser. It was painful. Horrifically painful.
The key failure mode was that, no matter how often I restarted/tuned, the etalon scan was always noisy.
The solution was to go into the user tuning mode for the Sprout and tune for power stability; this took ages (automated though)
Right afterwards, etalon looked nice and stable, and a quick scan worked nicely.

On the plus side, got a lot of coding done. The analysis script is now pretty plot-verbose and starting to look sensible.
Still no sign of transition?!
Maybe time to try another method; 
		* Tune up fibre coupling
		* Pull beam back
				* mount camera and observe spot size
		* Maybe before that we should switch back to BEC
				* Switch to TO with laser near transition
				* if no signal, tune closer to TO until signal acquired, then walk frequency back



## Hello!

* I ran a wide, fine run last night. Should be a ton of data. The analysis code is in good shape, too. Hopefully you can copy the data into a new dir, update the config with the new directory, and be good to go!
* If we don't see anything, I think we'll need a new method. We can't use a BEC as such, as it's not in the right state when trapped in the B trap.
* However, we CAN check alignment as we walk the laser away & kill the power, basically forming a dipole trap. Once we're close, we can switch to MOT and scratch our head about why we see nothing.
* Only after that should we try pulling a mirror back
* Checking out 0300 again, so excuse me if I'm a tad late tomorrow. What you can do:
		* Run the analysis; crack champagne?
		* If we're stashing the bubbly, walk back to tuneout (or just on this side of it) and get a trap freq signal back
	* Or, formulate another plan?!
* We're under a bit of a crunch now, we should put TO on the shelf until this laser goes home. 

## 2019-03-22, A New Hope
The overnight run died at some point because it lost Liquid N so need to add some hacks into the code to make the analysis work.

Salveging the TO sequnces 

		Intragation sequence : c:\remote\settings201906Mar124058.xml
		Calibration sequnce: c:\remote\settings201906Mar123950.xml

		c:\remote\settings201922Mar152500.xml


IT opt

init freq 
8.47		7.87		6.2
init amp
7.75		7			4.5
init_num
6439

Updating
freq 1
8.45		2158
8.46		4443
8.48		10602
8.49		13699	
8.50		13754
8.51		10041
8.52		4285
8.49		13231
amp 1
7.7			14030
7.60		14869
7.50		11890
7.40		10943
7.55		13063
freq 2
7.85		6036
7.86		6695
7.87		12908
7.88		14034
7.89		11473
a2
6.80		9866
6.90		6968
7.00		5929
7.10		9381
7.20		9537
6.70		6134

freq 1
8.52		6500
8.49		3500
8.51		6300
8.515		7200
freq 2
7.915		7200
7.89		6000
7.94		5800
7.96		2500
7.92		6500
7.90		5000

## (2019-03-25) Forward to BEC
New optimised TO settings
calibration : c:\remote\settings201925Mar141613.xml
interegation : c:\remote\settings201925Mar132223.xml

interegation trans: c:\remote\settings201925Mar143430.xml
c:\remote\settings201925Mar144631.xml
c:\remote\settings201925Mar151007.xml
c:\remote\settings201925Mar151532.xml

New Settings for transition:
calibration_file = 'c:\remote\settings201925Mar141613.xml';
interrogation_file = 'c:\remote\settings201925Mar151532.xml'; %1.4 V
interrogation_file = 'c:\remote\settings201925Mar155504.xml'; %0.7 V

## 2019-03-26

Transition sighted yesterday!


Things to do:
* Better logging
	* figure exporting
	* Generate metadata logfile 
* Desaturate signal
		* Lower optical power
				* Adjust set point
				* Expand beam
		* Higher init atom number
				* Visible with no evaporator?
* Narrower scan
* Calibration model
* Fitting
* First estimate!
* generalize monitor run feature


Comparing loading params

BEC sequence
LMOT frequency	
1		8.53
2		7.916
3		6.2
atten
17.8
2		7.1
3		4.5
Shunt
1		10
2		5.5
3		0
Quad
1		5.25
2		7.2
3		5.74

2sec transfer
quad
1		5.25		
2		7.2
Shunt
1		10
2		5.5

settings201925Mar141613 TO
settings201925Mar141613 BEC 
settings201927Mar180449 evap off, synth off
settings201927Mar184632 10 sec shorter				WORKS
settings201927Mar184948 AL attn fixed
settings201927Mar190600	

Moving to single-stage magnetic trap?
settings201927Mar191132.xml