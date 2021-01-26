## Standardised Set of Earthquake Records

There are three subfolders, t_0_02s, f_0_01s, t_0_005s, corresponding to three different time steps of 0.02s, 0.01s, and 0.005s respectively.

All records run for 100 seconds. The longest actual records lasts 80 seconds and all records are padded with zeros to ensure any transient behaviour is captured.
The shortsest records last only approximately 15 seconds, so contain significant zero-padding for the standard length of 100 seconds.

All of the .mat files contain two variables:

Time_vector  -  a row vector containing the time data for the chosen time step
Accel_matix  -  a matrix where each COLUMN is one earthquake record - either numbered 1 to 20 within each suite, or 1 to 60 for all records

All acceleration records are in m/s^2.

The Records data for each earthquake is:


    EQ# (overall) 		Suite		EQ# (within suite)		Location/Name					PGA (m/s^2)
	
     1			Medium		     1				Imperial Valley, (la01)    			4.5203
     2			Medium		     2				Imperial Valley, (la02)     			6.6288
     3			Medium		     3 				Imperial Valley, 1979, Array 5 (la03)     	3.8604
     4			Medium		     4				Imperial Valley, 1979, Array 5 (la04)     	4.7865
     5			Medium		     5				Imperial Valley, 1979, Array 6 (la05)     	2.9569
     6			Medium		     6 				Imperial Valley, 1979, Array 6 (la06)     	2.3008
     7			Medium		     7				Landers Eqk, 1992 (la07)     			4.1298
     8			Medium		     8				Landers Eqk, 1992 (la08)     			4.1749
     9			Medium		     9				Landers Eqk, 1992 (la09)     			5.097
    10			Medium		    10				Landers Eqk, 1992 (la10)     			3.5335
    11			Medium		    11				Loma Prieta, 1989, Gilroy (la11)     		6.5249
    12			Medium		    12				Loma Prieta, 1989, Gilroy (la12)     		9.5093
    13			Medium		    13				Northridge, 1994 (la13)     			6.6493
    14			Medium		    14				Northridge, 1994 (la14)     			6.4449
    15			Medium		    15				Northridge, 1994 (la15)     			5.233
    16			Medium		    16				Northridge, 1994 (la16)     			5.6858
    17			Medium		    17				Northridge, 1994, Sylmar (la17)     		5.5843
    18 			Medium		    18				Northridge, 1994, Sylmar (la18)     		8.0144
    19			Medium		    19				North Palm Springs, 1986 (la19)     		9.9943
    20			Medium 		    20				North Palm Springs, 1986 (la20)     		9.6761
------------------------------------------------------------------------------------------------------------------------------
    21			High 		     1				1995 Kobe (la21)     				12.58
    22			High 		     2				1995 Kobe (la22)     				9.0275
    23			High 		     3				1989 Loma Prieta (la23)     			4.0995
    24			High 		     4				1989 Loma Prieta (la24)     			4.6376
    25			High 		     5				1994 Northridge (la25)     			8.5162
    26			High 		     6				1994 Northridge (la26)     			9.2529
    27			High 		     7				1994 Northridge (la27)     			9.087
    28			High 		     8				1994 Northridge (la28)     			13.041
    29			High		     9				1974 Tabas (la29)     				7.9345
    30			High		    10				1974 Tabas (la30)     				9.7258
    31			High		    11				Elysian Park - Simulated (la31)     		12.712
    32			High		    12				Elysian Park - Simulated (la32)     		11.635
    33			High		    13				Elysian Park - Simulated (la33)     		7.6726
    34			High		    14				Elysian Park - Simulated (la34)     		6.6759
    35			High		    15				Elysian Park - Simulated (la35)     		9.7316
    36			High		    16				Elysian Park - Simulated (la36)     		10.793
    37			High		    17				Elysian Park - Simulated (la37)     		6.9784
    38			High		    18				Elysian Park - Simulated (la38)     		7.6131
    39			High		    19				Elysian Park - Simulated (la39)     		4.9058
    40			High		    20				Elysian Park - Simulated (la40)     		6.1328
------------------------------------------------------------------------------------------------------------------------------
    41			Low		     1				Coyote Lake, 1979, Gilroy Array 2 (la41)     	5.7834
    42			Low		     2				Coyote Lake, 1979, Gilroy Array 2 (la42)     	3.2681
    43			Low		     3 				Imperial Valley, El Centro Array 6 (la43)     	1.4067
    44			Low		     4				Imperial Valley, El Centro Array 6 (la44)     	1.0945
    45			Low		     5				Kern County Eqk, 1952 (la45)     		1.4149
    46			Low		     6				Kern County Eqk, 1952 (la46)     		1.5602
    47			Low		     7				Landers Eqk, 1992 (la47)     			3.3122
    48			Low		     8				Landers Eqk, 1992 (la48)     			3.0174
    49			Low		     9				Morgan Hill, 1984 Gilroy 3 (la49)     		3.1241
    50			Low		    10				Morgan Hill, 1984 Gilroy 3 (la50)     		5.3588
    51			Low		    11				Parkfield, 1966, Array 5 (la51)     		7.6565
    52			Low		    12				Parkfield, 1966, Array 5 (la52)     		6.1936
    53			Low		    13				Parkfield, 1966, Array 8 (la53)     		6.8001
    54			Low		    14				Parkfield, 1966, Array 8 (la54)     		7.7505
    55			Low		    15				North Palm Springs, 1986 (la55)     		5.0758
    56			Low		    16				North Palm Springs, 1986 (la56)     		3.7166
    57			Low		    17				San Fernando, 1971 (la57)     			2.4814
    58			Low		    18				San Fernando, 1971 (la58)     			2.2654
    59			Low		    19				Whittier, 1987 (la59)     			7.537
    60			Low		    20				Whittier, 1987 (la60)     			4.6907
