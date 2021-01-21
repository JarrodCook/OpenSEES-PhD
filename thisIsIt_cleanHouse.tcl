# THIS IS IT

# clear previous variables and commands
wipe

# establish model domain with 2 dimensions
model BasicBuilder -ndm 2

#**********Input parameters*****************************************

set LCfactor 20;

#*****Fixed or derived:*****

set pi [expr 2.0*asin(1.0)];	# Definition of pi
# puts "pi = $pi"
set g 9.81;						# gravity acc.

# Frame parameters
set Kground 1e11;				# ground stiffness
set Khorz 1e10;					# horizontal support stiffness
set A 0.1;						# area of frame elements
set E 200e9;
set M 1e5;						# frame mass
set W [expr $M*$g];				# weight force
set damping 0.03;				# critical damping ratio
set B 2.5;						# half width
set nodeMassx 1e3;				# nominal node mass - x dof
set nodeMassy 1e3;				# nominal node mass - y dof
set nodeMassz 0;				# nominal node mass - z dof

set MAll [expr ($LCfactor + 1) * $M];
set H [expr $aspectRat*2*$B];	# height
set HM [expr $H*2/3];			# mass height
set wn [expr 2*$pi/$Tn];		# natural frequency (fixed base)

set Kinitial [expr $wn*$wn*$MAll];	# initial system stiffness
set Icol [expr $Kinitial*$HM*$HM*$HM/(3*$E)];
# puts "Icol: [format {%0.2f} $Icol]"
set H2 [expr $H - $HM];
set Icol2 [expr $Kinitial*$H2*$H2*$H2/(3*$E)];

set IbeamNom [expr $Kinitial*$B*$B*$B/(3*$E)];		# inertia of beam elements
set Ibeam [expr 128*$IbeamNom]; # the beams have greater flexural stiffness than the columns, creates reasonable response
# puts "Ibeam: $Ibeam"
set Ibeam 10.0;

# Moment at uplift
set Mover [expr $Gup*$g*$MAll*$HM/$R];
set MW [expr $W*$B];
set MED [expr $Mover*(100 - $PTper)/100];
set MPT [expr $Mover - $MW - $MED];
set Fyielddev [expr $MED/(2*$B)];
set initialFPT [expr $MPT/$B];
set FyoverW [expr $Fyielddev/$W];

puts "Mover: [format {%0.4E} $Mover], MW: [format {%0.4E} $MW], MPT: [format {%0.4E} $MPT], MED: [format {%0.4E} $MED]"
# puts "MW: [format {%0.2E} $MW]"
# puts "MPT: [format {%0.2E} $MPT]"
# puts "MED: [format {%0.2E} $MED]"
puts "Fyield: [format {%0.4E} $Fyielddev], Fyield/W: [format {%0.2f} $FyoverW], PT force: [format {%0.4E} $initialFPT]"
# puts "Fyield/W: [format {%0.2f} $FyoverW]"
# puts "PT force: [format {%0.2E} $initialFPT]"


# Leaning column parameters
set L [expr 3*$B];				# distance to leaning column
set massLC [expr $LCfactor*$M];	# mass of leaning column
set WLC [expr $massLC*$g];		# weight force of leaning column
set WAll [expr $MAll*$g];		# weight force of full system

# Rigid link to LC
set Arigid 1000.0;				# define area of truss section (make much larger than A of frame elements)
set Irigid 100000.0;			# moment of inertia for p-delta columns  (make much larger than I of frame elements)
set Erigid 200e9;				# steel Young's modulus

# # Device and PT stiffness for uplift stiffness ratio alpha
# set alpha 0.8;
# set const1 [expr $HM/5/$B/$B];
# set EPTdev [expr $const1 * (4*$Kinitial*$HM + $WAll)];

set bGNG 2e-4;					# GNG post-yield stiffness ratio

# Device and PT stiffnesses for uplift stiffness ratios alpha and gamma
set alpha 0.8;
set gamma 0.1;
set const1 [expr $HM*$HM/$B/$B/4/(1-$bGNG)];
set const2 [expr ($alpha*$Kinitial)/(1-$alpha) - ($gamma*$Kinitial)/(1-$gamma)];
set EGNG [expr $const1*$const2];
set const3 [expr $HM/$B/$B];
set const4 [expr ($gamma*$Kinitial*$HM)/(1-$gamma) + $WAll];
set EPTstar [expr $const3*$const4 - 4*$bGNG*$EGNG];

# Device parameters
# set Edev [expr 20*$Kinitial];	# initial elastic modulus
# set Edev 2e8;
# set Edev $EPTdev;
set Edev $EGNG;


# PT parameters
# set EPT [expr 20*$Kinitial];	# elastic modulus
# set EPT 5e7;
# set EPT $EPTdev;
set EPT $EPTstar;
# no implied relation between Fy and W! (just that devices are at edges and PT is at centre)
# set initialFPT [expr ($Fyielddev*2*$PTper/(100 - $PTper)) - $WAll];
set initialStrainPT [expr -$initialFPT / $EPT];	# initial strain

# set Mbaserock [expr $W*$B + $initialFPT*$B + $Fyielddev*2*$B]
# puts "Mbaserock: [format {%0.2E} $Mbaserock]"

#**********Node placement*****************************************

# Rocking frame
#		no.	x	y
node	1	0	0;		# base centre (frame)
node	2	0	$HM;	# mass height (frame)
node	3	-$B	0;		# left foot (frame)
node	4	$B	0;		# right foot (frame)
node	5	-$B	0;		# left foot (ground)
node	6	$B	0;		# right foot (ground)
node	7	0	0;		# base centre (ground)
node	10	0	$H;		# roof (frame)

# Leaning column
#		no.	x	y
node	8	$L	0;		# base of leaning column
node	9	$L	$HM;	# roof of leaning column

#**********Nodal masses*****************************************

# Rocking frame
#		node	x			y			z
mass	1		$nodeMassx	$nodeMassy	$nodeMassz;		# bottom centre (frame)
mass	2		$M			$M			$nodeMassz;		# mass height
mass	3		$nodeMassx	$nodeMassy	$nodeMassz;		# left foot
mass	4		$nodeMassx	$nodeMassy	$nodeMassz;		# right foot
mass	10		$nodeMassx	$nodeMassy	$nodeMassz;		# roof

# Leaning column
#		node	x			y			z
mass	8		0			0			$nodeMassz;		# base of leaning column (pinned)
mass	9		$massLC		$massLC		$nodeMassz;		# top of leaning column (mass)

#**********Fixity conditions*****************************************

# Rocking frame
#	node	x	y	z
fix	5		1	1	1;	# ground below left foot - fully fixed
fix	6 		1	1	1;	# ground below right foot - fully fixed
fix	7 		1	1	1;	# ground below bottom centre - fully fixed

# Leaning column
#	node	x	y	z
fix 8 		1	1	0;	# base of leaning column - pinned

#**********Geometric tranfromations*****************************************

# Tags
set transfTag_C 1; # column tag
set transfTag_B 2; # beam tag

# PDelta
geomTransf PDelta $transfTag_C;	#0 0 -1
geomTransf PDelta $transfTag_B;	#0 1 0

# Linear
# geomTransf Linear $transfTag_C;
# geomTransf Linear $transfTag_B;

#**********Material models*****************************************

# high compressive stiffness rocking edge
set rockMatTag 1
uniaxialMaterial ENT $rockMatTag $Kground

# elastic post-tensioning
set PTMatTag 2
uniaxialMaterial Elastic $PTMatTag $EPT

# high compressive stiffness horizontal support at rocking edges
set supportMatTag 3
uniaxialMaterial ENT $supportMatTag $Khorz

# nominal small stiffness elastic material to provide stability in zero stiffness plastic cases
set nominalKMatTag 25
uniaxialMaterial Elastic $nominalKMatTag 1

# device hysteresis GNG - LHS
set GNGLHSMatTag 4
uniaxialMaterial GNG $GNGLHSMatTag $Edev $Fyielddev $Pitch $bGNG
set GNGLHSrockMatTag 5
uniaxialMaterial Parallel $GNGLHSrockMatTag $rockMatTag $GNGLHSMatTag $nominalKMatTag

# device hysteresis GNG - RHS
set GNGRHSMatTag 6
uniaxialMaterial GNG $GNGRHSMatTag $Edev $Fyielddev $Pitch $bGNG
set GNGRHSrockMatTag 7
uniaxialMaterial Parallel $GNGRHSrockMatTag $rockMatTag $GNGRHSMatTag $nominalKMatTag

# elastic post-tensioning with inital strain (large yield strains to avoid plastic behaviour)
set PTisMatTag 8
uniaxialMaterial ElasticPP $PTisMatTag $EPT 1000 -1000 $initialStrainPT

# define truss material for link to leaning column
set TrussMatTag 9
uniaxialMaterial Elastic $TrussMatTag $Erigid;

#**********Elements*****************************************

# Frame elements
element elasticBeamColumn 1		1 2		$A $E $Icol		$transfTag_C;	# vertical frame element - below mass
element elasticBeamColumn 11	2 10 	$A $E $Icol2	$transfTag_C;	# vertical frame element - above mass
element elasticBeamColumn 2		1 3 	$A $E $Ibeam	$transfTag_B;	# left horizontal frame element
element elasticBeamColumn 3		1 4 	$A $E $Ibeam	$transfTag_B;	# right horizontal frame element

# element elasticBeamColumn 2		1 3 	$Arigid $Erigid $Irigid	$transfTag_B;	# left horizontal frame element
# element elasticBeamColumn 3		1 4 	$Arigid $Erigid $Irigid	$transfTag_B;	# right horizontal frame element

# Rocking edge elements
element zeroLength 4 5 3 -mat $GNGLHSrockMatTag -dir 2;	# left rocking edge (vertical)
element zeroLength 5 6 4 -mat $GNGRHSrockMatTag -dir 2;	# right rocking edge (vertical)
element zeroLength 6 5 3 -mat $supportMatTag	-dir 1;	# left rocking edge (horizontal)
element zeroLength 7 6 4 -mat $supportMatTag	-dir 1 -orient -1 0 0 0 -1 0;	# right rocking edge (horizontal)

# Base centre element
element zeroLength 8 7 1 -mat $PTisMatTag 		-dir 2;		# PT element

# Rigid link to leaning column
element truss  9 2 9 $Arigid $TrussMatTag;	# rigid pinned element between the 2 masses (produces cleaner curve than equalDOF command)
# equalDOF 2 9 1; # node slaved in x-direction to roof node of frame 

# Leaning column
element elasticBeamColumn  10  8  9 $Arigid $Erigid $Irigid $transfTag_C;	# rigid pinned column

#**********Gravity loading*****************************************

# Gravity loading
pattern Plain 1 Linear {
    load 2 0 -$W 0;		# weight force applied to top point of frame (-y)
	load 9 0 -$WLC 0;	# weight force applied to leaning column
}

#**********Analysis commands*****************************************

constraints Plain;				# constraint equation setting
numberer Plain;					# numbering scheme used to assemble the system of equations
system BandGeneral;				# system of equations
algorithm Linear;				# solution algorithm
integrator LoadControl 0.1;		# incremental solution via load control (steps of 0.1 * load pattern)
analysis Static;				# type of analysis
set ok [analyze 10];			# analyze 10 steps (of 0.1 * load pattern)
loadConst -time 0;				# sets loads constant and resets times to 0.0 (gravity load always applied)
puts "Gravity load applied"

#**********Eigenvalues*****************************************

set Eigenvalue [eigen -fullGenLapack 1]
puts "Eigenvalue : [format {%0.4E} $Eigenvalue]"

#**********Damping properties*****************************************

set omega1 [expr {sqrt(abs($Eigenvalue))}];	# high frequency bound (root of first eigenvalue - fundamental period)
set omega2 [expr {$omega1/10}];	# low frequency bound (10 * fundamental period)
set alphaM [expr {2*$damping*$omega1*$omega2/($omega1 + $omega2)}];	# mass proportional damping term
set betaK 0;			# stiffness proportional damping term
set betaKinit 0;		# initial stiffness proportional damping term
set betaKcomm [expr {2*$damping/($omega1 + $omega2)}];	# committed stiffness proportional damping term
rayleigh $alphaM $betaK $betaKinit $betaKcomm ;# Rayleigh damping command

#********************************Outputs**********************************

set KrigidrockWiebe [expr $EPT*$B*$B - $WAll*$HM]
# puts "Krigidrock: $Krigidrock"
set KrockWiebe [expr 1 / (($HM*$HM/$KrigidrockWiebe) + (1/$Kinitial))]
puts "KrockWiebe = [format {%0.4E} $KrockWiebe]"

set Krigidrock [expr $EPT*$B*$B + $Edev*(2*$B)*(2*$B) - $WAll*$HM]
# puts "Krigidrock: $Krigidrock"
set Krock [expr 1 / (($HM*$HM/$Krigidrock) + (1/$Kinitial))]

set KrigidrockYield [expr $EPT*$B*$B + $bGNG*$Edev*(2*$B)*(2*$B) - $WAll*$HM]
# puts "KrigidrockYield: $KrigidrockYield"
set KrockYield [expr 1 / (($HM*$HM/$KrigidrockYield) + (1/$Kinitial))]

puts "Kinitial: [format {%0.4E} $Kinitial], Krock: [format {%0.4E} $Krock], KrockYield: [format {%0.4E} $KrockYield]"
# puts "Krock: [format {%0.2E} $Krock]"
# puts "KrockYield: [format {%0.2E} $KrockYield]"

# set Kdev [expr ($HM*$HM*$Kinitial + $WAll*$HM) / ($B*$B + 2*$B*2*$B)]
# puts "Kdev: [format {%0.4E} $Kdev]"

set xuplift [expr $Mover/$HM/$Kinitial]
puts "Approximate mass deflection at uplift: $xuplift"

set xyield [expr $xuplift + $aspectRat*$Fyielddev/$Edev*2/3]
puts "Approximate mass deflection at yield: $xyield"

