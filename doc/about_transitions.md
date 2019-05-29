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
	* Performed some open-heart surgery on the cavity while measuring power at fibre input. Looks like beam was clipping on the case lid on the way out! Shuffling it a couple mm provided a beautiful beam spot and a nice high power throughput. Some alignment later and we have fibre coupling of 7.5/18 = 42# at low power, up to 80/220 (36#) at high power

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
interrogation_file = 'c:\remote\settings201925Mar151532.xml'; #1.4 V
interrogation_file = 'c:\remote\settings201925Mar155504.xml'; #0.7 V

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


## 2019-03-28
Some sequences built while looking for differences between 0.6 and 2s transfer
I am suspicious that the drifts affecting the optical molasses stage could account for our loss of fluro over time, molasses stage has not been attended to in some time. This is not likely to be the cause of our lost BEC, though, that was probably our fault.

15 sec b trap test settings201927Mar224917	
extended 0.6s hold settings201928Mar005945 
2sf pared back to minimal load settings201928Mar011156 

## 2019-03-29

30K ATOMS IN BEC YAY YAY YAY

Laser locked, BEC looking good

Setting up for search and scan
		* Peak previously sighted at 175MHz high from 727303244.58893MHz,  the Ritz predicted frequency
		* That is, peak at 727303414.58893 MHz 
		* We pass through an AO which passes the +1 order of a 189 MHz modulation, so to compensate for this we must offset by 189/2 MHz in the red
		* Therefore this peak should be centred at **363651612.794465** MHz

calibration settings201929Mar012806
measurement settings201929Mar012745
Probe power:
		140mW into AOM
		74mW at fibre coupler
		26mW after first filter cube :(

Laser reaches 200mV set point on photodiode for 19.5mW after waveplates. Neat conversion factor, hey.

Ok. Let's look for signal, probe beam set point 0.15V

Sighting at 		363.651693
Saturated by 		363.651730

BEC vanishes and remaining atoms really hot! Perfect

Bingo, saturates at 363651708 MHz
Remove support post from lens translation mount; still saturated!
Move final telescope lens all the way away from chamber; STILL saturated
Cut set point to .05V settings201929Mar024626, still saturated
Take final lens OUT, beam now GIGANTIC.
SIGNAL NOT SATURATED. Barely. Well, still does, but not every time.
Scan manually; 
Saturates between 363651707 and 363651699
-> Centre at 363651703
Shorten exposure by 2 seconds; only on during MOT load
Still saturates at 363651702.5 MHz
Hmmm. Could be WM drift at this stage. Looks pretty stable shot to shot, but can't say without calibration
FINALLY get a handle on saturation maybe? Unless it's drifted
Now have exposure only for 200ms during first ITC, with .05V set point and giant beam with photodiode gain up 10dB!
That's a factor of 20 in duration, 3 in amplitude, 10 in gain, plus an absolutely colossal beam to get away from saturation.


## Good morning!

Things we learned
		* This thing is super duper sensitive, we'll have no trouble saturating. Which might make it hard to control!
		* Also, we're applying the light in the ITC; I tried applying to MOT, much weaker effect. 
		* Beam appears to be desaturated and transition narrowed down to a few MHz! Would be good to run the code over the scan I left running.
Things to do:
		* Bring up atom number in BEC while running analysis on the morning run
		* Clean one desk or optics bench
		* Work towards computing the Allan deviation of the measured transition peak. Later we'll extend this to characterize the wavemeter with respect to the Cs transition.
		* When I'm in, let's calibrate the wavemeter and try get something running over the weekend. 
		* I'm going away this weekend and will need to prep once I've got some sleep, so I won't be in until after lunch.
Things to think about:
		* There's an ambient magnetic field that will Zeeman shift the transition. Could estimate this effect, I anticipate only a few MHz.
		* Is it worth concocting another sequence to obtain the zero-field value?


* Raised atom count by about 5k, not great but not bad
* Doing some comparison scans for the over night run, to see if it drifts/ is effected by other sources such as the evap
* added some features to the analysis code see side computer, value for the transition is very close to theory
* Should try to understand zeeman shift/ characterise it then move on to the next transition

## 2019-04-02
Calibrating wavemeter

Stop WM feedback					DONE
note setpoint 						372.196857/744.393723, tisaf 805.2
tisaf setpt to 822.2 				DONE		
run wm feedback without blue 		DONE	
run measure_2photon_transitions 	DONE
turn on pmt 						DONE
run measure 3 times 				DONE 	had to increase scan range, drift was approx 10MHz
stash in folder 					DONE
stop WM feedback 					DONE
turn on modulation fn 				DONE
select ch1 on settings signal 		DONE?
stop wm monitor						DONE
Operation->Calibration->etc 		DONE?
kill matlab   						DONE
Start wm monitor 					DONE
run ws8 fb then meas_2p    			DONE
Kill meas_2p 						DONE
reboot ws8 feedback  			 	DONE	
restore noted setpt 				DONE
power off pmt and mod 				DONE

Pre cal offset -6.125 MHz
Post cal offset 0.25MHz


## 2019-04-03

Issues calibrating wavemeter were thanks to incomplete to-do list!
Including the missing stage:

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
Run lock_to_2p
select ch1 on settings signal
stop wm monitor
Operation>Calibration>etc0
kill matlab 
run ws8 fb then meas_2p 
Kill meas_2p, reboot ws8 feedback & restore noted setpt
power off pmt & mod

So, after calibrating wavemeter again, we have offset calculations:
0.263787±0.029643
0.274663±0.033188
0.383103±0.037893
0.405642±0.026936
0.396748±0.041161
0.393202±0.034740
MEAN  0.35285758972168 STD 0.0652732911713087
Could probably fine tune a bit, but we're not running overnight so I'll leave laser locked overnight and see what happens come morning.


##2019-04-04

ZEEMAN DAY (or two or...)

* Run trap simulator with ITC settings stage 1 and 2

* Measure ITC AOM drive in first stages of ITC

* Use Mathematica script to estimate splitting in 5^3S_1 manifold transitions, use simulator and/or AO offset to estimate field 

* Write script (mathematica or matlab) to compute polz state (relative to lab) given waveplate angles and translate into B-field frame

	* The dream: A mathematica script that accepts a set of WP angles plus a magnetic field and predicts the line positions and relative strengths ;)

* Measure some lines in the 5^3S_1 manifold in first stage ITC

* If field still nice and homogeneous in stage 2 based on trap simulator, use second stage as high-field measurement

* 


Objective: One Good Run
Target is the 412nm transition to 5^3 S_1

Calibrated wavemeter 10:00
Measured offset from 3 scans of F3->F3 (check which) Cs transition:
Mean 0.373769693076611 STD 0.0705489367563666

Retuned laser to 824nm, recoupled fibre
180mW into AOM (could probably improve)
73mW into fibre
29mW after first post-fibre beam cube

Lock offset: 252.988722MHz when attached to PC

TODO: 
	Approximate beam profile (eg aperture over power meter, measure peak power and estimate 1/e^2 radius)
	Calc expected Zeeman shifts -> MOT AO, trap sim
	Find transition again
	Search for Zeeman transitions
	Test polarization dependence
	Try high-field measurement

peak around (record in the red):
first: 363651706.794468 (HWP = 146)
Switch to polarization of opposite handedness and search around the opposite side of the vacuum transition; nothing shows up.
An absence of extra peaks is consistent with the pump beam spin-polarizing the ground state to the m_j = 2 state. This can only scatter by sigma- transition to the m_j = 1 state of the 5^3 S_1; so we may actually be unable to resolve more Zeeman peaks! 
* We can make the basic check by scanning very wide around the sighted transition; if no lines are found, this supports the hypothesis that we are spin-polarized.
* We could try to pump into more states, but require a model of how the pump beam polarizes the atoms.
* We could try looking at another transition with more observable splittings; for example driving with linear probe beam on the 5^3D_5 transition we could resolve two peaks.
* We can try to characterize the action of the waveplates based on our measurements of their fast axes and test consistency with this hypothesis.


### Evening session
	* OK, let's pinpoint the peak and ensure consistency between control and analysis.
	* Run one scan 20MHz, 1MHz steps either side of the 5^3S_1 peak with HWP at 150 degrees (very circular)
		Fitted single peak at 727303224.46MHz
		2x Scan centre - 189MHz is 727303224.588936!!! Looks like we've got agreement at last
	* Run another scan with waveplate at 190 degrees - should be linear
		* Only one peak sighted, within 0.67MHz error of previous run. Expected amplitude to decrease, but hard to say given 1MHz steps easily missed the peak value.


Got locked out of lab, progress significantly stifled :(

* Setting up for long scan
	* If this is the 2^3P_2 m_j=2 to 5^3S_1 m_j=1, then it is expected to be approx 20MHz red of the vacuum resonance. The greatest splitting is approx 74MHz at 18G between the delta_j = +1 and -1 transitions from the g_ground = 0 state. For safety, I set up a scan 100MHz wide centred approx on the vacuum wavelength, 727303244 (which actually is just about bang-on the theoretical value)
	* I set the qwp to 190 degrees to pass linear (aka ambidextrous) light into the chamber, so if we are able to drive other sigma transitions we push them both
	* A quick test shows that we can still see the peak of the transition even though half the light is polarized in the other handedness, so we'll see at least one.
	* Another quick test: Set the probe to switch on during second stage cooling, settings201904Apr231624
		* Expect peak to move. If AO calibration is to be believed, it should shift 9MHz to the blue, closer to theory.
		* Fitted peak at 727303233 MHz, shifted blue by 9MHz YES!!! Physics works!
	* Probably nothing else to learn near this transition, so...


* Tune laser to 805.4nm
* 15 mW delivered after fibre beam cube
* 47 mW input to fibre
* 91 mW input to AOM
* Laser output power 132 mW :( after optimizing the crap out of the cavities - Try different crystal position if we need more power
* Note that while walking laser, kept eye on post-fibre photodiode and tweaked alignment on the run
* Predict peaks at 744396232.78, 744396178.2, and 744396186.92 MHz should be visible*
* so set scan centred at 744396205MHz (midpoint of extremal values each 27MHz away
	* Set scan range 40MHz either side (in the blue), so 80MHz wide in the blue in 1MHz steps
	* Each scan should take 2.3hrs so we can run 3-4 before analysing. While we wait, plenty to work on.
* Set exposure time to 300ms to make up for loss of power and splitting polz between sigma+ and sigma-


## 2019-04-05
Today's tasks:
Leave the scan running as long as you can bear. You should be able to work out a method for task 2 given data in this log/wiki, and it should be pretty straightforward.
	* First: Read wiki. Measure WM offset. Use the closest reference available. More than one if practical. Then clean/tidy something.
	* Central goal for today: Back-calculate field as precisely as possible, with quantified uncertainty. 
	* Extremely useful but not urgent: complete reference wavelength table with closest Cs references to all observed/target wavelengths, list red and blue freqs, sorted for readability
	* Bonus: build PD set point check into analysis program
	* Bonus: Add intensity calibration to analysis script to provide vertical scale
In future:
	* Chase down branching ratios/Einstein A coefficients from forbidden transition
	* Write a function that turns waveplate dial angles into polarization in the lab/atom frame (also useful for TO, we can fix the handedness here)

Running analysis on over night run

average offset 0.263610861518166 std 0.278591590217379

## Weekend scan

LET'S GET CRAZY
Biggest distance between transitions to the triplet D state is 332MHz (744396179.59MHz to 744396512.46 MHz) So I will scan 200MHz (blue) either side of the midpoint 744396346MHz. That's 5.5 hours per scan with 1MHz steps. This should do the trick!
I'll run the interrogation period in the second stage of ITC so we get a second field value for the data collected today. 
* PD set point 350mV, 14.4mW, QWP 146 degrees, fresh dewar in

Things to do:
* Saturday:		Come in, change waveplate angle, repeat scan
* Sunday: 		Come in, run scan in first stage of ITC with your choice of waveplate angle
* Monday: 		Party?

* Saturday: 
	* Left beam dump in, so no data :( But at least we have big, stable BECs. 
	* Measuring WM offset:
		mean -0.006408±0.000345 (se), sd 0.004145
	* Take the block out, restart run.
	* Came back at end of run to check output & change parameters
		* Landscape looks pretty barren. Would probably not lose much by restricting scan to near peaks. However, for completeness, let's get 'em all for now.
		* Modified interface so it cycles through multiple probe settings; in this case the sequence is calibrate-> Stage 1 -> stage 2
		* This way I can run the scan overnight and get both magnetic field settings for a single WP angle, will take 8.2hrs per scan
		* Set QWP to 236, so we'll get
			* Stage 1 D_2, qwp 236 <- already have this
			* Stage 1 D_3, qwp 236 <- 
			* Stage 2 D_2, qwp 236  |-- New data
			* Stage 2 D_3, qwp 236 <-
		* Note that an extra segment is required in the analysis script to deal with these 
	* So why do we only see these peaks? I wonder; can we actually see the D_3 and the D_1, because the D_2 is forbidden by selection rules?! 
	* Some of these transitions are very close to theory. We NEED to calibrate wavemeter in a closer region.
	* Also, what about peaks that are too weak for us to recover?
	* NB QWP fast axis offset at 203 degrees, HWP fast axis offset at 20 degrees
* Sunday
	* Damnit, unexpected failure; Huge numbers of coincidences in LabView log times?! How did that happen? Probably salvageable 
	* FUCKING DAYLIGHT SAVINGS. We should probably query a server somewhere for UNIX time. Or disable auto-DST updates?
	* Hack fix: Mask out all the files in that 2hr window...
	* Setting up for another run. What do we need?
	* Confident that the shift in peak amplitude is enough to identify everything. So take waveplates out and run two-stage scan
	* Goal: Determine magnetic field. Along the way, will likely need to identify transitions. So, here goes.
		* Step 1: Extend analysis capability DONE
		* Step 2: Identify lines...
			* Almost there!

## Hello, Kieran!

Obviously, quite a lot happened on the weekend, so I'm going to take a chunk of Monday to stay home. Some progress to acquaint yourself with:
	* the zeeman_theory code in the repo is pretty well established (could be more flexible), have a play with it. Uses include:
		* Determining scan ranges
		* Generating theory to fit to data
		* Getting a feel for the Zeeman effect
		* The plots are pretty poorly labeled atm, but if the colour scheming works then:
			* darker = lower initial ground m_j, lighter = higher ground m_j 
			* ^ = sigma+, o = pi, v = sigma-
	* The plot in the lower right of main monitor is about as far as I got. It's almost done! Try manually identifying the observed peaks with the predicted data, then fit the theory offset & magnetic field strengths.
		* Note the displayed data came from me ripping the waveplates out! 
		* We also have a ton of data with various light polarizations
		* Sorry about the poor code hygiene/half-baked figures...
	* Nonetheless, I think it's pretty concievable that you can make serious headway towards these today:
		* What is the magnetic field strength?
		* What pumped states (in the 2^3P manifold) exist in the trap?
		* Which transitions are we driving? 
	* I'll come in for a bit in the afternoon, but you can contact me on messenger any time.
	* Other sundry things that need attention:
		* Fill up Little Boy, Fat Man has been running all weekend
		* Find Cs transitions that we could use to calibrate near 402 (and the forbidden line, fingers crossed)
		* Fix the caching functionality in the analysis script; often re-imports when not necessary, taking heaps of time
		* Add the analog import back in, and add a check for the PD voltage
			* If the PD voltage is below set point but steady, we can just scale the signal at that shot, if the signal is linear...
		* So let's find out! Take a few quick scans over a peak with varying intensity and see how the peak height varies.Probably the best data to take now, and quite quick too. 
		* Code could do with a general tidy & better documentation
		* Good luck and happy science!

## 2019-04-09
Measured offset of  -3.00047796169917 p\m 0.8 MHz AM
					-0.823342±0.001718 PM

Setting up WM stability analysis;
LOCKED to f_cs_2p_6SF3_8SF3 = 364507238.363 MHz
	Calibrate WM 364.5072384THZ @ 1940
	Left to run overnight.
	IN MORNING
		* Stash WM logs
		* Measure 2p transition f_cs_2p_6SF3_8SF3, AND at least two others, record offsets.
		* Calibrate if drift is real bad
		* Go about your day :)

## 2019-04-10
Overnight wm loggs stashed
Measured offset f_cs_2p_6SF3_8SF3: 1.181198±0.001235
Measured offset f_cs_2p_6SF4_8SF4: 2.330434±0.005532

c:\remote\settings201910Apr100304.xml

Stage 1: c:\remote\settings201910Apr100911.xml (old c:\remote\settings201909Apr182710.xml)
c:\remote\settings201910Apr123035.xml
TO DO:
5_3D1 signal stregth with qwp angle
Table of B fields for all stages
Cs Calibration around 804 nm
Prove there's only one peak for 5_3D1
Good scan for the 5_3D1 peak


Updated in-trap cooling settings settings201910Apr160535
### Setting up for RF spectroscopy 
	* Targets: Measure the minimal trap splitting at equivalent settings to the first and second stages of in-trap cooling:
		* Stage 1: Shunt 10V,  quad 5.25V
		* Stage 2: Shunt 5.5V, quad 7.2V
	* Starting from BEC settings, pretty cold but only ~10k, settings201910Apr160535, start building slow ramp up to final field value
		* settings201910Apr161812 quad field ramps to 5.25V in 1sec after evap, AL attn maxed during pulse
			* Nothing seen; RF shield left on during ramp, probably killing BEC 
		* Ok, nothing reappears when RF shield disabled. May as well work towards target incrementally then, if sudden change doesn't work.
		* settings201910Apr162901 1sec ramp to quad=3.6V
		* 
	* NB: PAL centred on 1.739MHz, and our estimates from trap sim/AOM/observed Zeeman shifts put the field values at around 10 and 14 Gauss (plus or minus two in each case), so we can begin our sweeps centred at 14 and 25MHz, respectively.
	* Atom number crashed pretty hard, had to realign some optics. Back to approx 45k atoms no sweat, super chill. TRY AGAIN TOMORROW
	* Tune laser to 3D1 transition on first stage, hits it right away. Let's see what we can do with the waveplate.
	* Signal saturates with waveplate at 146 degrees, vanishes at 236 degrees!
	* Therefore we ascertain the handedness; 236 degrees is sigma+, 146 is sigma-. Signal barely under saturation when set to midpoint 190 degrees.
	* Let's double-check our initial state overnight. Set the waveplate to 190 degrees, so hopefully a gross mixture of light. 
		* Oh, laser was set by Tenma the whole time, so massively overexposed and broadened for sure. So, return to LV set point.
		* Signal at about half full dynamic range when at 190 degrees, perfect
	* Predicted lines from the 2^3P_2_mj_1 state appear between 744.396440 THz and 744.396500 THz, giving a 60MHz scan range in the blue either side of the mean 744.396470. Including the m_g=0 state we predict transitions up to 744.396530 THz, so I set up a scan 45MHz either side of 744.396485 MHz. In 1MHz steps (blue), that's 90 steps per scan, 2 measurement & 1 cal = 1.8 hours per scan. 
	* If BEC number hold up, we should have a good half dozen scans by the morning which should let us pin down even little peaks
	* WM not calibrated because this isn't really a 'measurement' run, just a hypothesis test.
	* Scan started at 1930


## 2019-04-11

Double trouble overnight: Seed laser offset drifted and killed BEC. Also, doubler had a hard time locking after about an hour.
	* Tune up etalon and cavity lock error signal after rebooting control unit (modulated noise on etalon lock error). 
	* Offset barely lasted an hour, goddamn. BEC comes back when offset retuned... 
	* Alright fine, let's run again, but stay to babysit. During runtime, we should fix up the code to:
		* Check whether blue laser is locked
			* Reject shots with too great a difference between the blue and 2x red, or if WM error too great
		* Check photodiode voltage & reject if too noisy. For the moment, we can log act/set value and then scale signal, until linearity verified
		* Reject if atom number too low
	* Misc code jobs:
		* Fix caching?!
		* Implement SLOGS
		* Break config into master/slave

 * CHECKLIST FOR GOODLY DATA RUNS
 	* Logging
 		* WM
 			* WM receiving blue & red light
 			* WM_main has useblue enabled
 		* Analog import
 			* Analog USB input is running
 			* AI trigger is 100ms before the probe AOM switches on
 			* Analog USB capture time 200ms longer than exposure time
		* Hardware setup
			* Probe output & lock signals optimized at centre of scan range
			* Probe aligned to reach at least 20# above desired setpoint
			* Measure probe power before second post-fibre beamcube 
 	* Calibrate wavemeter & measure nearby transitions
 	* Restart DLD acquistion & clear DLD output director
y 	* hit auto-run & free run
 	* Babysit to ensure:
 		* Analog import updating
 		* Probe sequence changing between cal/stage 1/stage 2
			* WM set point updating
			* Atom number stable
	* Create local_config.m with fields
		* pd_setpoint
		* pd_delay = 100
		* wp_angle
		* transition_level
		* Other custom settings

## 2019-04-12

RF spectroscopy sequence c:\remote\settings201912Apr121027.xml

c:\remote\settings201912Apr135548.xml

c:\remote\settings201912Apr144045.xml

c:\remote\settings201912Apr145344.xml

next version using the synth c:\remote\settings201912Apr160418.xml


## 2019-04-15 WM error curve
cs_2p_6SF3_8SF3: Running mean 1.529590±0.007065 (se), sd 0.035323
cs_2p_6SF4_8SF4: Running mean 2.834866±0.000963 (se), sd 0.004814
cs_6SF3_6PF2co3: very large (can't observe)

## 2019-04-18

Short day, starts at 1500.
Goals for today & tomorrow
* Polish code to presentation plots
* Look for calibration lines
* Analyse linearity data
* Fine up ITC spectroscopy
* Final scans of 5^3D and 5^1S
* Search for new transitions overnight
* Compute expected gain for forbidden transition
* Lock to Cs and record WM stability over the weekend

Goals for after weekend:
* Go straight for Forbidden transition, we're running out of time

* LabView being strange... Uh-oh.
In the meantime, I'll try this calibration business. Don't need LV for that.
Pre-calibration measurements

transition 			Offset
cs_2p_6SF4_8SF4 	300kHz 
CALIBRATED

Looking for others:
transition 			Offset
cs_2p_6SF3_8SF3 	-1.35 MHz
cs_2p_6SF4_8SF4 	-54kHz
Can't find others, and adding all these switch cases is unwieldy. I'll make the Cs transitions into an easily navigable struct, then call it a day and crush out some goodness tomorrow. 
Left laser locked to cs_2p_6SF4_8SF4, will examine log in the morning (assuming it stays locked oops)

## 2019-04-19

WM offset on the F4-F4 transition is  1.03 MHz
WM offset on the F3-F3 transition is -0.41 MHz
Try a few others, can't spot 'em...
	* Power too low?
	* Cell alignment lost when retuning etalon?
Perhaps calibrate and leave measure running over weekend to look at long-term wm drift
Moving on to measuring stuff...
BEC happy :) 50k counts
Look for 5^3D1 and 5^3S today

Scanning over 5^3D1
Stage 2 power: 0.15V/5.19mW after PD cube, beam fully defocused (use previous estimate of radius)
Stage 1 power: 0.1V/3.4mW
qwp 146 degrees
Data collected, 667 shots
Argh, was only running stage 1.
Set up again to scan over both, then will leave overnight. Will likely die over the weekend. RIP


## 2019-04-23

Setting up for overnight run of the 5^3D transitions
Run over the weekend had setpoints up in the 700THz regime so unlocked early on and died :(
Kieran tweaked up atom number and recorded WM offsets
Jacob takes over
	Calibrate wavemeter on F4-F4
	Can't measure F3_F3, laser keeps unlocking
	Restart control unit after finding that etalon lock error offset incorrect
	Tune laser back to operating region
	At 0.1V set point, power delivered is 4.85mW

Tasks for Wednesday:
	* Look for Cesium reference lines
		* Turn on PMT and roll through Cs transitions to spot as many as we can find. 
	* Check independence from probe beam parameters
		* Scan reference peak (say, 5^3S1?)
		* Repeat for pump power & detuning
			* Adjust ITC power/freq as much as possible while retaining BEC
				* Measure IT cooling power if practical. If not, measure pre-ND filter and use that to normalize
			* Measure difference in power/freq
			* Re-scan peak
	* Set up for overnight scan
	* 

## 2019-04-24
The run taken after resetting the offset last night shows some nice looking peaks around the D_2, D_3 region, but very small. Goes to show their relative strengths are quite different. For the same intensity, the D_1 is barely visible, so perhaps we'll have to run them separately after all.
Some issues with the code to fix:
	* Number culling DONE
	* Probe setpoint check
		* Serves as proxy for ECD lock check also
	* Implemented, but it's pretty rough atm
		There is an AC background (5kHz?) which has not been subtracted
		Amplitude is approx 5mV so artificially increases PD range by about 10mV
		Data is smoothed almost immediately, so leads to some error in calculating on-time
		Offset on photodiode signal implies that power is not linear wrt. absolute voltage
	* Probably sufficient for now, let's move on to other tasks
Calibration curve revisited
	Measuring WM offsets:
		F4-F4 -0.696162±0.006263 (se), sd 0.031313
		F3-F3 -2.519053±0.026838 (se), sd 0.134191
	Calibrate to F3
		F4-F4 1.907405±0.002462 (se), sd 0.012308
		F3-F3 0.410344±0.001532 (se), sd 0.007659 <- expected, error less than precision of calibration


Sight no more. Return more deliberately another time.

FINALLY get a good scan of the 5^3D_1, with tightly focused lens and a mere 10ms exposure time. 
Signal is about half of full dynamic range, ideal to explore a few dependencies. The transition is pretty narrow so scans should be pretty quick.
PLAN
* Take fast scan around 5^3D_1, stage 1 ITC
* Adjust ITC detuning and re-scan
* Adjust ITC power and re-scan
* Adjust probe power and re-scan
* Retake original scan to check drift

Hm, can't find spectrum analyser!!
However, do have measurements of AO offset from magnetic field calibrations:
	138MHz and 133.6MHz during first stages of ITC, which have voltages approx 8.5 and 7.9. 
		Therefore estimate dependence of approx 7.6MHz/V.
			Annoyingly, we can't really move the ITC this much! 
	Oops, initial scan had steps too large to be practical.
	So I'll come in tomorrow (ANZAC) and finish this up.
	Meantime, may as well take the real good scan.

## 2019-04-25
Made a few improvements to code
	analog checks now implemented
	Local configurations now set in local_opts.m file in data directory
	Zeeman correction implemented kinda roughly
		Should functionalize and accept a string labeling excited state, then calculating automatically
	Changed signal to atom number loss, but pretty easy to reconfigure to use ratio (see calibration section)
	Text output now a bit more complete
	output saves automatically and more figures saved

As at 2141:
F3-F3 offset -0.030973±0.005948 (se), sd 0.029742
F4-F4 offset 1.519872±0.007246 (se), sd 0.036231
Reinstalled waveplate to allow WM to record blue WL & maintain blue lock
10.8mW post-PD cube
Large beam
Testing linearity/probe stark
250ms exposure time, 4cm beam waist
set point 	Post-PD power 	peak centre 			peak height 	peak ratio 		Peak width (MHz)
0.10V 		3.31mW			744396515.790(0.264) MHz    12508 			0.36			1.67(0.48)
0.15V 		5.14mW 			744396516.282(0.180) MHz 	19248			0.62			2.20(0.41)
0.20v 		6.5mW  			744396516.508(0.128) MHz	18604 			0.74			2.28(0.16)
0.25V 		8.9mW  			744396515.106(0.157) MHz	23202			0.84			2.37(0.35)
Recorded PD ranges: [0.0475,0.06450.0817,0.099] resp, with '0' reading 0.011

All peaks heights are differences in atom number
Peak ratios are normalized: 1-(calibration-N_atom_probe_on)/calibration
All within wavemeter accuracy, and I've seen bigger centre variations between successive Cs measurements. 
Conclusion: No observable shifts. Running at high power is good.
	Notice that the WM lock is a bit noisier? Even when running with the blue lock disabled. 
	1: Need to stabilize, perhaps reboot laser
	2: Take a more detailed look at the WM set determination; perhaps average over the period where the probe is high

OK, what else? 
Need good run of 5^3D_1, 5^3D_2,3, 5^3S_1 (?), look for 4^1D_2 (4^1P_1?)
- The run we now have of the 5^1D_2 with 1MHZ steps looks great; Zeeman corrected values agree to 0.7MHz
So I'm happy to call that done.
- So, try to get a nice picture of the D_2 and D_3, then tomorrow:

## 2019-04-26
Tasks:
	* Ensure that overnight run went OK
	* Try to change pump detuning by at least a MHz, and see whether this shift a peak. That is;
		* Pick a nice peak (5^3D_1?) and take a pretty coarse scan. I found that scan_range = 6 and stepsize 0.5MHz (red) is pretty quick and provides a good fit, and it's easy to operate near saturation.
		* Kick the ITC settings around. You can try to calculate the expected shift in probe freq to get started, but eventually it must be measured with the spectrum analyser. 
			* The ITC intensity might be too hard to measure directly, but you could measure the power before the ND filter and just note the proportional power throughput. Try measuring at the trap anyhow, that'd be nice to know.
	* In the end, we have to try these, but they might be impractical, which we'll have to live with.
	* Fine scan of 5^3S_1
	* Look for 4^1* lines/set up scan
	* Figure out why WM lock is so unstable now?! There's a worrying amount of drift in between updates
		* Can  remedy by more detailed checking of WM set point in code, but should fix the problem at the root too!
Goal: These transitions DONE by monday. Next week, forbidden search. 


Oh yeah, Bryce and student are doing ZS velocity distribution measurements today, so not a lot of opportunity to take data.

Last night's run worked OK, but I forgot to include the second stage probe, oops... 
Good news is that the side peaks are JUST visible, and could fit to them, but they're real small. Nice and distinct too. I think we can take another scan with more power, presuming the atom number can be bolstered a bit (just saturated overnight), or just live with a saturated main peak to bring out the rest. 


More detailed WM log retrieval:
AI log written at end of recording period
LV log written at start of run
So, probe_on_time = ai_time - (ai_duration - pd_delay)
all computable from ai log, then passes a posix interval to wm log which is sampled & averaged

## 2019-04-28
	Leaving temp control unit unplugged overnight seems to fix the problem, laser now nice and stable again
	Running with large beam, 0.25V set point, 8.45mW after PD cube, 60k atoms yea yea
	Cs transitions;
	cs_2p_6SF3_8SF3 	-3.041977±0.004047 (se), sd 0.008095
	cs_2p_6SF4_8SF4 	-1.554317±0.001144 (se), sd 0.003433
	Calibrated to F4
	cs_2p_6SF4_8SF4 	0.062275±0.002331 (se), sd 0.011653
	cs_2p_6SF3_8SF3 	-1.342083±0.000165 (se), sd 0.000827
	Would be interesting to observe effect of calibration on transition peak centres....


* Oops, saturated. Well, at least we'll get a look at where the peaks are. Can re-scan if really necessary. 
* I ran at 250mV, which I reckon might not have saturated had the number stayed steady. Oh well.

Some things to do;
	* Run local_opts in master_config instead of adding stacks of files to the path
	* Pepper some comments through code
	* Fix caching...
	* Cull points that are too far outside LV set bounds, presumably in bin_by_wavelength
		* WTF is with those errand points if the bins are too small?


	* Transition final checklist:
		name 		Directory
		5^3D_1 		20190425_5^3D1_overnight (No analog logs... but beautiful!)
		5^3D_2,3	20190429_cal_5^3D_2_3_qwp_142
		5^1D_2 		20190417_5^1D2_mj_1_itc_both_qwp_146_overnight (also no analog logs, but very pretty...)
		5^3S_1 		20190429_5^3S_1_qwp_146_both_stage   		7.2mW
		4^1D_2 		
		4^3D_?

	Plus have a lot of prelim stuff (without full logging) to bolster stat precision, but there's not much to win there
	Tune laser back to 824nm, obtain 14.3mW power post-PD at 680mV (20dB gain, has been for ages)



## 2019-05-09

### Stat errors   
Level				mean MHz 					theory diff MHz 	Past value (them-us) THz
5^1D2 	 			744430345.471(0.047)		2.121				N/A
5^3S1 				727303247.812(0.143)		3.212				727.31696(.013713)
5^3D1  				744396515.588(0.224)    	4.448				744.41009(.013575)
5^3D3 	            744396204.115(0.341)       -4.245               744.41009(.013886)
					744396202.6225(0.6) 	   -5.737				Including all 4 peaks
5^3D2 				744396235.844(0.282) 		8.264 				744.41009(.013855)

Note that past expt could not distinguish between the lines between P2,1 and D1,2,3!
Nor the P2,1 and 3S1 lines

Error includes fit error and Zeeman error
Does not include setpt err (laser noise)

### Sys error budget

Probe power 	2MHz ultra conservative
WM drift 		
WM accuracy 	20MHz :(



Uncertainty in freq set point?
Quick notes on final measurements
name 	 power 	 		last WM cal offset (F3, F4) 			 theory difference		
5^1D2 	 				1.529590±0.007065,2.834866±0.000963      
											300kHz next day
5^3S1 	 7.2mW?			-1.342083±0.000165,0.062275±0.002331
5^3D1 	                -0.030973±0.005948, 1.519872±0.007246
5^3D2,3	 8.45mW			-1.342083±0.000165,0.062275±0.002331

Powers for first two: Reload LV settings and check set points, then use power measurements to ballpark
Beam radius for 5^1D is small

#### Probe power dep
From analysis (see power_dependence.m), dependence apparently swamped by drift
WORST CASE analysis: assume WM steady and fit trend to lower points => 0.225MHz/mW
.225MHz/mW, or 852 Hz/(mW/cm^2). 
Should remake plot with new ZS analysis
Calculate in terms of intensity and estimate worst case for 5^1D
#### Zeeman error
Can go back and find splittings more accurately
In the meantime, error in freq(1MHz) => fractional error in B => fractional error in calculated Z shift 


Systematics:
	Pump effects => Too delicate to measure, disadvantage of technique
	WM error 	 => last calibration + drift analysis
			 	 => Ultimately, manufacturer spec... **20MHz** in blue!


How to interpred the D2,3 data?

Peaks stage 1:
[744396159.014 744396190.235 744396221.123]
stage 2:
[744396180.504 744396204.973 744396222.922]

Hm, puzzling when compared to theory. Looks like two sigma- transitions. But the other; pi or sigma+?

The WP was actually set to 142 degrees, but that's still not much sigma+. 
Ok, well, if the transitions are strong enough, we can't have seen sigma+ anywhere else
And there's no trace of pi in the other observation (though it was weak af); did we scan far enough?
At least in some, yeah... Soo what if we leave that peak alone and just treat the ones we understand?
Target states in ascending freq order:
5^3D_3_2, 5^3D_2_1, 5^3D_3_3?!

Welp. Maybe even the middle peak is hokum, but we can get something for the 3D_3. The others, well, we can use the 3D_2 to fix the expected 

Pump beam detunings:
g_g = 2
g_e = 1.5
d_mg = 1
mu = 1.4MHz/G

So in stage 1, cooling transition is split 25.55 MHz, and 16.002 in stage 2
Therefore the beams are (ish) detuned 0.5 and 2MHz respectively

Ok. So peak figures are pretty well at acceptable draft level, but need to be combined into one for display, probably. 
Next thing to do is write a post-process function that imports results and fig data, then combines them nicely.

Other figure to-do:
	* Level diagram
	* Expt setup
