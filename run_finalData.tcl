# THE FINAL ANALYSIS

#**********File settings*****************************************
#
# Save file settings
set recorderdir "finalResults";# recorder directory name
file mkdir $recorderdir ;# create results directory
set GMdir "GM/AllEQs/"

set f [open scaleFactors1.txt r];
set scaleFactors1 [split [string trim [read $f]]];
close $f;

# Force at uplift in g (multiplied by Wtrib by H (and divided by R) to get moment at uplift)
set f [open Gups.txt r];
set Gups [split [string trim [read $f]]];
close $f;

set maxPeriods 100;
set periodStep 0.05;

# Create result files
set testResults [open $recorderdir/testResults.txt w]
# Node displacements
set nodeD [open $recorderdir/nodeD.txt w]
# Node reactions
set nodeR [open $recorderdir/nodeR.txt w]
# Node accelerations
set nodeA [open $recorderdir/nodeA.txt w]

#**********Time settings*****************************************

set dt_analysis 1e-3;	# time step for transient analysis
set max_time 100.0;		# maximum time

#**********Loading pattern*****************************************

# Reduction factor
set R 8;

# Set uplift PT contribution (percentage)
set PTper 66;

# Select ground motion records
set GMFiles [glob -nocomplain -directory $GMdir -type f *.txt];

# Define aspect ratios
set aspectRats [list 2 4 6 8];
set NARs [llength $aspectRats];

# Define natural periods (seconds)
set Tns {{0.2 0.3 0.4} {0.4 0.5 0.6 0.7} {0.5 0.6 0.7 0.8 0.9 1.0} {0.6 0.7 0.8 0.9 1.0 1.1 1.2}};
set NTns [llength $Tns];

# Define Pitches (metres)
set Pitchs [list 1e-3 2e-3 5e-3 10e-3 20e-3];

# Preallocate analyses counter
set theCtr 0;

# Start the clock
set startTAll [clock seconds]

# Loop ground motion records
foreach gMotion $GMFiles {
	
	# Set current ground motion record
	set gMotionName [string range $gMotion 0 end-4 ]
	set gMotionNumber [string range $gMotion 12 end-4 ]
	
	# Loop structure aspect ratios
	for {set i 0} {$i < $NARs} {incr i 1} {
		
		# Set current aspect ratio
		set aspectRat [lindex $aspectRats $i]
		
		foreach Tn [lindex $Tns $i] {
		
			# Set uplift force
			set indexRef [expr round($Tn/$periodStep) -1]; # index for scale factor and uplift force lists
			set Gup [lindex $Gups $indexRef]; # Uplift force in g
			
			# Set scale factor
			set ScaleFactorRef [expr (($gMotionNumber-1)*$maxPeriods) + $indexRef -1];
			set scaleFactor [lindex $scaleFactors1 $ScaleFactorRef];
			
			# Loop device pitches
			foreach Pitch $Pitchs {
				
				# Load the model and run gravity analysis
				source finalData.tcl
				
				# ---earthquake ground motion---
				timeSeries Path 2 -dt $dt_analysis -filePath $gMotion -factor $scaleFactor -prependZero
				pattern UniformExcitation 2 1 -accel 2
				
				#**********Recorders**********
				
				# Node displacements
				recorder EnvelopeNode -file $recorderdir/nodeDisps.txt	-node 1 2 3 4 9 10 -dof 1 2 3 disp
				# Node reaction forces
				recorder EnvelopeNode -file $recorderdir/nodeReactions.txt	-node 1 2 3 4 5 6 7 8 9 10 -dof 1 2 reaction
				# Node accelerations
				recorder EnvelopeNode -file $recorderdir/nodeAccels.txt	-node 1 2 3 4 9 10 -dof 1 2 accel
				
				#**********Analysis commands*****************************************

				constraints Plain;					# constraint equation setting
				numberer RCM;						# numbering scheme used to assemble the system of equations
				system BandGeneral;					# system of equations
				test EnergyIncr 1e-12 100 0;		# convergence test
				algorithm Newton;					# solution algorithm
				integrator Newmark 0.5 0.25;		# Newmark Beta integration scheme with constant average acceleration
				analysis Transient;					# type of analysis

				# run the analysis
				set startT [clock seconds]
				set ok [analyze [expr {int($max_time/$dt_analysis/2)}] $dt_analysis]
				puts "50% complete"
				set ok [analyze [expr {int($max_time/$dt_analysis/2)}] $dt_analysis]
				puts "Done!"
				
				set endT [clock seconds]
				puts "Execution time: [expr $endT-$startT] seconds."

				#**********Display commands*****************************************

				# Write results to file
				set demandLeft [eleResponse 4 material 1 material 2 demand]
				set demandRight [eleResponse 5 material 1 material 2 demand]
				set ratLeft [eleResponse 4 material 1 material 2 ratchetCount]
				set ratRight [eleResponse 5 material 1 material 2 ratchetCount]
				puts $testResults "$gMotionNumber $Tn $scaleFactor $aspectRat $Pitch $PTper [format {%0.4E} $Fyielddev] $ok [format {%0.2f} $Eigenvalue] $demandLeft $demandRight $ratLeft $ratRight"
				
				wipe;
				
				#**********Write analysis recorder results to analyses results files
				
				# Node displacements
				set f [open $recorderdir/nodeDisps.txt];
				set nodeDNow [split [string trim [read $f]] "\n"];
				close $f;
				puts $nodeD "[lindex $nodeDNow end]"
				
				# Node reaction forces
				set f [open $recorderdir/nodeReactions.txt];
				set nodeRNow [split [string trim [read $f]] "\n"];
				close $f;
				puts $nodeR "[lindex $nodeRNow end]"
				
				# Node accelerations
				set f [open $recorderdir/nodeAccels.txt];
				set nodeANow [split [string trim [read $f]] "\n"];
				close $f;
				puts $nodeA "[lindex $nodeANow end]"
				
				# Update counter
				set theCtr [expr $theCtr+1];
			}
		}
	}
}

set endTAll [clock seconds]
puts "Total execution time for $theCtr analyses: [expr $endTAll-$startTAll] seconds."

close $testResults
close $nodeD
close $nodeR
close $nodeA

wipe ;# this will kill the recorders which will cause them to close the files

