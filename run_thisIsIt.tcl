# THIS IS IT


#**********File settings*****************************************

# Save file settings
set recorderdir Analysis_Results/thisIsIt ;# recorder directory name
file mkdir $recorderdir ;# create results directory
set GMfile "earthquake signals/Standard Set/t_0_001s/EQ1.txt"

set EQNumber [string range $GMfile 43 end-4 ]
# puts "EQ number: $EQNumber"

set f [open scaleFactors1.txt r];
set scaleFactors1 [split [string trim [read $f]]];
close $f;

set maxPeriods 100;
set dTn 0.05;

# Force at uplift in g (multiplied by Wtrib by H (and divided by R) to get moment at uplift)
set f [open Gups.txt r];
set Gups [split [string trim [read $f]]];
close $f;

# Save analysis details file
set testDetails [open $recorderdir/testDetails.txt w]

#**********Time settings*****************************************

set dt_analysis 1e-3 ;# time step for transient analysis
set max_time 100.0 ;# maximum time

#**********Model settings*****************************************

# Reduction factor
set R 4;

# Frame
set Tn 0.7;			# fixed base period
set aspectRat 4;	# aspect ratio

# Devices
set Pitch 5e-3;		# tooth pich
set PTper 100;		# % of force from PT at uplift

# Derived
set indexRef [expr round($Tn/$dTn) -1];
# puts "indexRef: $indexRef"
set Gup [lindex $Gups $indexRef]; # Uplift force in g
# puts "Gup: $Gup"

set ScaleFactorRef [expr (($EQNumber-1)*$maxPeriods) + $indexRef -1];
set scaleFactor [lindex $scaleFactors1 $ScaleFactorRef];
# puts "scale factor: $scaleFactor"

# set scaleFactor 1.0; #**************************************************************************

# Load model file
# source thisIsIt_cleanHouse.tcl
source thisIsIt_cleanHouseNoDevices.tcl

#**********Loading pattern*****************************************

# ---initial velocity loading (via empty UniformExcitation)---
# timeSeries Constant 3 -factor 0.0
# pattern UniformExcitation 2 1 -accel 3 -vel0 7.0

# ---earthquake ground motion---
timeSeries Path 2 -dt $dt_analysis -filePath $GMfile -factor $scaleFactor -prependZero
pattern UniformExcitation 2 1 -accel 2

# ---acceleration ramp---
# timeSeries Path 3 -dt $dt_analysis -filePath AccelRamp.txt -factor 1.0 -prependZero
# pattern UniformExcitation 3 1 -accel 3

#**********Recorders*****************************************

# node displacements
recorder Node -file $recorderdir/Disps.txt -time -node 1 2 3 4 5 6 7 8 9 10 -dof 1 2 3 disp
# node accelerations
recorder Node -file $recorderdir/Accs.txt -node 1 2 3 4 5 6 7 8 9 10 -dof 1 2 3 accel
# node reactions
recorder Node -file $recorderdir/Reactions.txt -node 1 2 3 4 5 6 7 8 9 10 -dof 1 2 3 reaction
# element forces
recorder Element -file $recorderdir/Forces.txt -ele 1 2 3 4 5 6 7 8 9 10 11 force
# element deformations
recorder Element -file $recorderdir/Defs.txt -ele 1 2 3 4 5 6 7 8 9 10 11 deformation
# demand
recorder Element -file $recorderdir/Demand.txt -ele 4 5 material 1 material 2 demand
# ratchet count
recorder Element -file $recorderdir/Ratchet.txt -ele 4 5 material 1 material 2 ratchetCount

recorder EnvelopeNode -file $recorderdir/envelopeDisps.txt	-node 1 2 3 4 9 10 -dof 1 disp

#**********Analysis commands*****************************************

constraints Plain;					# constraint equation setting
# constraints Penalty 1.e15 1.e15
numberer RCM;						# numbering scheme used to assemble the system of equations
system BandGeneral;					# system of equations
test EnergyIncr 1e-12 100 0;		# convergence test
# test NormDispIncr 1.0e-6 10 0
algorithm Newton;					# solution algorithm
integrator Newmark 0.5 0.25;		# Newmark Beta integration scheme with constant average acceleration
analysis Transient;					# type of analysis

# run the analysis
set startT [clock seconds]
set ok [analyze [expr {int($max_time/$dt_analysis/4)}] $dt_analysis]
puts "25% complete"
set ok [analyze [expr {int($max_time/$dt_analysis/4)}] $dt_analysis]
puts "50% complete"
set ok [analyze [expr {int($max_time/$dt_analysis/4)}] $dt_analysis]
puts "75% complete"
set ok [analyze [expr {int($max_time/$dt_analysis/4)}] $dt_analysis]
puts "Done!"


puts "OK: $ok"	;# will output 0 if it worked, other number if did not


set endT [clock seconds]
puts "Execution time: [expr $endT-$startT] seconds."

# write results to file
puts $testDetails "$EQNumber $Tn $scaleFactor $Gup $aspectRat $Pitch $PTper $R [format {%0.4E} $Fyielddev] $ok [format {%0.2f} $Eigenvalue]"

wipe;

# read peak displacements
set f [open $recorderdir/envelopeDisps.txt];
set maxDisps [split [string trim [read $f]] "\n"];
close $f;
puts "Peak displacements, at mass: [lindex $maxDisps 1 1], at top: [lindex $maxDisps 1 5]"

# pattern Plain 2 Linear {
    # load 2 1 0 0
# }

# set ctrl_node 2
# set ctrl_dof 1
# set max_drift 0.10
# set max_disp [expr {$max_drift*$H}]
# set nsteps 10000
# set disp_incr [expr {$max_disp/$nsteps}]

# test EnergyIncr 1e-12 100 0
# algorithm Newton
# integrator DisplacementControl $ctrl_node $ctrl_dof $disp_incr
# analysis Static

# set startT [clock seconds]
# set ok [analyze $nsteps]
# set endT [clock seconds]

#**********Display commands*****************************************

close $testDetails

wipe ;# this will kill the recorders which will cause them to close the files

