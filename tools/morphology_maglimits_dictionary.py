
def maglimits(x):
    x = str(x)
    out = {
        '5122' : [15., 20.],
        '5123' : [15., 20.],
        '5124' : [15., 20.6],
        '5125' : [15., 20.],
        '5126' : [15., 19.6],
        '5127' : [15., 20.],
        '5128' : [15., 20.],
        '5129' : [15., 20.],
        '5081' : [15., 20.],
        '10774' : [16., 20.],
        '3102' : [15., 20.],
        '3103' : [15., 20.],
        '8540' : [15., 20.],
        '8750' : [15., 20.],
        '8751' : [15., 20.],
        '8752' : [15., 20.],
        '8753' : [15., 20.],
        '3637' : [15., 20.],
        '3638' : [15., 20.],
        '3131' : [15., 20.],
        '3132' : [15., 20.],
        '8765' : [15., 20.],
        '2628' : [15., 20.],
        '2629' : [15., 20.],
        '9571' : [15., 20.],
        '9572' : [15., 19.3],
        '5214' : [17., 19.8],
        '5215' : [15., 20.],
        '2672' : [15., 20.],
        '9576' : [15., 19.6],
        '8818' : [15., 19.3],
        '8309' : [15., 20.],
        '8310' : [15., 20.],
        '3191' : [15., 20.],
        '3192' : [15., 20.],
        '8313' : [15., 19.1],
        '8314' : [15., 20.],
        '9578' : [15., 20.],
        '8384' : [15., 20.],
        '9858' : [15., 19.6],
        '9579' : [15., 20.],
        '5254' : [17., 19.6],
        '9353' : [15., 18.3],
        '3726' : [15., 20.],
        '3727' : [15., 20.],
        '3728' : [15., 20.],
        '3729' : [15., 20.],
        '9582' : [15., 19.3],
        '1691' : [15., 20.],
        '9884' : [15., 19.3],
        '4765' : [15., 20.],
        '4769' : [15., 20.],
        '9378' : [15., 19.3],
        '9379' : [15., 20.],
        '9380' : [15., 20.],
        '9381' : [15., 20.],
        '9382' : [15., 20.],
        '4776' : [15., 20.],
        '4782' : [15., 20.],
        '10013' : [16., 20.],
        '4784' : [15., 20.],
        '3253' : [15., 20.],
        '3254' : [15., 20.],
        '4792' : [15., 20.],
        '9588' : [15., 20.],
        '10015' : [16., 20.],
        '8381' : [15., 20.],
        '4798' : [15., 20.],
        '8383' : [15., 19.1],
        '4800' : [15., 20.],
        '10102' : [16.5, 20.],
        '3270' : [15., 20.],
        '10103' : [16., 20.],
        '3278' : [15., 20.],
        '5069' : [15., 19.],
        '8312' : [15., 20.],
        '8483' : [15., 20.],
        '3288' : [15., 20.],
        '3805' : [15., 20.],
        '3806' : [15., 20.],
        '3807' : [15., 20.],
        '9954' : [15., 19.1],
        '9955' : [15., 20.],
        '4842' : [15., 20.],
        '4843' : [15., 20.],
        '4844' : [15., 20.],
        '4845' : [15., 20.],
        '4846' : [15., 20.],
        '4847' : [15., 19.5],
        '4848' : [15., 20.],
        '1265' : [15., 20.],
        '3314' : [15., 20.],
        '8436' : [15., 20.],
        '9978' : [15., 20.],
        '9979' : [15., 20.],
        '8444' : [15., 20.],
        '3331' : [15., 20.],
        '3334' : [15., 20.],
        '8461' : [15., 20.],
        '2833' : [15., 20.],
        '3346' : [15., 20.],
        '9859' : [15., 19.3],
        '3350' : [15., 20.],
        '3351' : [15., 20.],
        '3352' : [15., 20.],
        '5082' : [15., 20.],
        '10014' : [16., 20.],
        '2847' : [15., 20.],
        '2848' : [15., 20.],
        '2849' : [15., 20.],
        '5083' : [15., 19.1],
        '9181' : [15., 20.],
        '8497' : [15., 20.],
        '9182' : [15., 20.],
        '8505' : [15., 20.],
        '9183' : [15., 20.],
        '9184' : [15., 20.],
        '8516' : [15., 20.],
        '10053' : [17., 20.1],
        '10054' : [17.5, 20.3],
        '9580' : [15., 20.],
        '2889' : [15., 20.],
        '2890' : [15., 20.],
        '9186' : [15., 20.],
        '8382' : [15., 20.],
        '8528' : [15., 20.],
        '10104' : [17.5, 20.3],
        '4948' : [15., 20.],
        '4750' : [15., 20.],
        '4950' : [15., 20.],
        '4951' : [15., 20.],
        '4952' : [15., 20.],
        '4953' : [15., 20.],
        '4954' : [15., 20.],
        '4955' : [15., 20.],
        '4956' : [15., 20.],
        '4957' : [15., 20.],
        '8472' : [15., 20.],
        '2915' : [15., 20.],
        '2916' : [15., 20.],
        '9573' : [15., 19.6],
        '9574' : [15., 20.],
        '8552' : [15., 20.],
        '9577' : [15., 19.6],
        '3434' : [15., 20.],
        '3435' : [15., 20.],
        '3436' : [15., 20.],
        '3437' : [15., 20.],
        '3438' : [15., 20.],
        '9583' : [15., 20.],
        '9584' : [15., 18.6],
        '9585' : [15., 19.6],
        '9586' : [15., 20.],
        '9587' : [15., 20.],
        '8564' : [15., 20.],
        '9589' : [15., 20.],
        '9590' : [15., 20.],
        '9591' : [15., 20.],
        '9592' : [15., 20.],
        '9593' : [15., 20.],
        '9594' : [15., 20.],
        '9918' : [15., 20.],
        '10089' : [17.5, 20.1],
        '1938' : [15., 19.6],
        '3475' : [15., 20.],
        '1954' : [15., 19.6],
        '9883' : [15., 20.],
        '3313' : [15., 20.],
        '10055' : [17., 20.1],
        '3229' : [15., 20.],
        '3507' : [15., 20.],
        '3508' : [15., 20.],
        '3509' : [15., 20.],
        '9156' : [15., 19.1],
        '9157' : [15., 20.],
        '2507' : [15., 20.],
        '2508' : [15., 20.],
        '2509' : [15., 20.],
        '2511' : [15., 20.],
        '2512' : [15., 20.],
        '2514' : [15., 20.],
        '2516' : [15., 20.],
        '2517' : [15., 20.],
        '2518' : [15., 20.],
        '2519' : [15., 20.],
        '2520' : [15., 20.],
        '2521' : [15., 18.6],
        '2522' : [15., 20.],
        '2523' : [15., 20.],
        '2524' : [15., 20.],
        '2525' : [15., 20.],
        '2526' : [15., 18.8],
        '2527' : [15., 20.],
        '2528' : [15., 19.6],
        '9185' : [15., 20.],
        '2530' : [15., 20.],
        '2531' : [15., 20.],
        '2535' : [15., 20.],
        '2536' : [15., 20.],
        '9980' : [15., 20.],
        '2033' : [15., 20.],
        '2035' : [15., 20.],
        '9581' : [15., 20.],
        '8311' : [15., 20.],
        '9575' : [15., 20.], #up to now UPAD-DATA-ACCESS and EDR tiles
        '9254' : [15., 19.1], #from now on UPAD TEST-2 tiles (T2)
        '9255' : [15., 20.],
        '9256' : [15., 19.3],
        '9257' : [15., 20.],
        '9258' : [15., 19.6],
        '9259' : [15., 18.8],
        '9260' : [15., 19.3],
        '9261' : [15., 20.],
        '8973' : [15., 20.],
        '8840' : [15., 20.],
        '8853' : [15., 20.],
        '8865' : [15., 20.],
        '5853' : [15., 20.],
        '8877' : [15., 20.],
        '8889' : [15., 20.],
        '5819' : [15., 20.],
        '5820' : [15., 20.],
        '5821' : [15., 20.],
        '5822' : [15., 20.],
        '5823' : [15., 20.],
        '5824' : [15., 20.],
        '5825' : [15., 20.],
        '5826' : [15., 20.],
        '8901' : [15., 20.],
        '5832' : [15., 20.],
        '5833' : [15., 20.],
        '5834' : [15., 20.],
        '5835' : [15., 20.],
        '5836' : [15., 20.],
        '5839' : [15., 20.],
        '5840' : [15., 20.],
        '8914' : [15., 20.],
        '5843' : [15., 20.],
        '5844' : [15., 20.],
        '5845' : [15., 20.],
        '5846' : [15., 20.],
        '5847' : [15., 20.],
        '5848' : [15., 20.],
        '5850' : [15., 20.],
        '5851' : [15., 20.],
        '5852' : [15., 19.6],
        '8925' : [15., 20.],
        '5854' : [15., 20.],
        '5855' : [15., 20.],
        '5856' : [15., 20.],
        '5857' : [15., 20.],
        '5858' : [15., 20.],
        '5859' : [15., 20.],
        '5861' : [15., 20.],
        '5862' : [15., 20.],
        '5863' : [15., 20.],
        '5864' : [15., 20.],
        '5865' : [15., 20.],
        '5866' : [15., 20.],
        '5867' : [15., 20.],
        '5868' : [15., 20.],
        '5869' : [15., 20.],
        '5870' : [17., 20.],
        '5871' : [15., 20.],
        '5872' : [16.6, 20.],
        '5873' : [15., 19.3],
        '5874' : [15., 20.],
        '5875' : [15., 20.],
        '5876' : [15., 20.],
        '5877' : [15., 20.],
        '5878' : [15., 20.],
        '5879' : [15., 20.],
        '5880' : [15., 20.],
        '5881' : [15., 20.],
        '5882' : [15., 20.],
        '5883' : [15., 20.],
        '5884' : [15., 20.],
        '5885' : [15., 20.],
        '5886' : [15., 20.],
        '5887' : [15., 20.],
        '5888' : [15., 20.],
        '8961' : [15., 20.],
        '5890' : [15., 20.],
        '5891' : [15., 20.],
        '5892' : [15., 20.],
        '5893' : [15., 20.],
        '5894' : [15., 18.8],
        '5895' : [15., 20.],
        '5898' : [15., 20.],
        '5899' : [15., 20.],
        '5900' : [15., 18.],
        '5901' : [15., 20.],
        '5902' : [15., 20.],
        '5903' : [15., 20.],
        '5904' : [15., 20.],
        '5905' : [15., 20.],
        '5906' : [15., 20.],
        '5907' : [15., 18.8],
        '5908' : [15., 20.],
        '5909' : [15., 20.],
        '5910' : [15., 20.],
        '5911' : [15., 20.],
        '5912' : [15., 20.],
        '5913' : [15., 20.],
        '5914' : [15., 20.],
        '5915' : [15., 20.],
        '8989' : [15., 20.],
        '9001' : [15., 20.],
        '5931' : [15., 20.],
        '5932' : [15., 20.],
        '5933' : [15., 20.],
        '5934' : [15., 20.],
        '5935' : [15., 20.],
        '5936' : [15., 20.],
        '5937' : [15., 20.],
        '5938' : [15., 20.],
        '5939' : [15., 20.],
        '5940' : [15., 20.],
        '5941' : [15., 20.],
        '5942' : [15., 20.],
        '5943' : [15., 20.],
        '5944' : [15., 20.],
        '5945' : [15., 20.],
        '5946' : [15., 20.],
        '5948' : [15., 20.],
        '5949' : [15., 20.],
        '5950' : [15., 20.],
        '5951' : [15., 20.],
        '5952' : [15., 20.],
        '9025' : [15., 19.3],
        '5954' : [15., 20.],
        '5955' : [15., 20.],
        '5956' : [15., 20.],
        '9014' : [15., 20.],
        '5958' : [15., 20.],
        '5959' : [15., 20.],
        '5960' : [15., 20.],
        '5961' : [15., 20.],
        '9037' : [15., 20.],
        '5842' : [15., 20.],
        '9049' : [15., 20.],
        '9060' : [15., 19.1],
        '9073' : [15., 20.],
        '8937' : [15., 20.],
        '9084' : [15., 19.3],
        '9365' : [15., 18.3],
        '5953' : [15., 20.],
        '9097' : [15., 19.3],
        '9109' : [16., 20.],
        '9121' : [15., 20.],
        '9132' : [15., 20.],
        '9145' : [15., 20.],
        '8949' : [15., 20.],
        '11264' : [15., 20.], #from now on UPAD TEST-3 tiles (T3)
        '11265' : [15., 20.],
        '11266' : [15., 20.],
        '11267' : [15., 20.],
        '11268' : [15., 20.],
        '11269' : [15., 20.],
        '11270' : [15., 20.],
        '11271' : [15., 20.],
        '11272' : [15., 20.],
        '11273' : [15., 20.],
        '11274' : [15., 20.],
        '11275' : [15., 20.],
        '11276' : [15., 20.],
        '11277' : [15., 20.],
        '11278' : [15., 20.],
        '11279' : [15., 20.],
        '11280' : [15., 20.],
        '11281' : [15., 20.],
        '11282' : [15., 20.],
        '11283' : [15., 20.],
        '11284' : [15., 20.],
        '11285' : [15., 20.],
        '11287' : [15., 20.],
        '11288' : [15., 20.],
        '11289' : [15., 20.],
        '11290' : [15., 20.],
        '11291' : [15., 20.],
        '11292' : [15., 20.],
        '11293' : [15., 20.],
        '11295' : [15., 20.],
        '11296' : [15., 20.],
        '11297' : [15., 20.],
        '11298' : [15., 20.],
        '11299' : [15., 20.],
        '11300' : [15., 20.],
        '11301' : [15., 20.],
        '11302' : [15., 20.],
        '11303' : [15., 20.],
        '11305' : [15., 20.],
        '11306' : [15., 20.],
        '11307' : [15., 20.],
        '11308' : [15., 20.],
        '11309' : [15., 20.],
        '11310' : [15., 20.],
        '11311' : [15., 20.],
        '11312' : [15., 20.],
        '11313' : [15., 20.],
        '11314' : [15., 20.],
        '11315' : [15., 20.],
        '11316' : [15., 20.],
        '11317' : [15., 20.],
        '11318' : [15., 20.],
        '11319' : [15., 20.],
        '11320' : [15., 20.],
        '11321' : [15., 20.],
        '11322' : [15., 20.],
        '11323' : [15., 20.],
        '11324' : [15., 20.],
        '11325' : [15., 20.],
        '11326' : [15., 20.],
        '11327' : [15., 20.],
        '11329' : [15., 20.],
        '11330' : [15., 20.],
        '11331' : [15., 20.],
        '11332' : [15., 20.],
        '11333' : [15., 20.],
        '11334' : [15., 20.],
        '11335' : [15., 20.],
        '11336' : [15., 20.],
        '11337' : [15., 20.],
        '11340' : [15., 20.],
        '11341' : [15., 20.],
        '11342' : [15., 20.],
        '11343' : [15., 20.],
        '11344' : [15., 20.],
        '11345' : [15., 20.],
        '11347' : [15., 20.],
        '11348' : [15., 20.],
        '11349' : [15., 20.],
        '11350' : [15., 18.],
        '11351' : [15., 20.],
        '11352' : [15., 20.],
        '11353' : [15., 20.],
        '11354' : [15., 20.],
        '11355' : [15., 20.],
        '11356' : [15., 20.],
        '11357' : [15., 20.],
        '11358' : [15., 20.],
        '11359' : [15., 20.],
        '11360' : [15., 20.],
        '11361' : [15., 20.],
        '11362' : [15., 20.],
        '11363' : [15., 20.],
        '11364' : [15., 20.],
        '11365' : [17., 19.],
        '11366' : [15., 20.],
        '11367' : [16.5, 18.6],
        '11368' : [15., 20.],
        '11369' : [15., 20.],
        '11370' : [15., 20.],
        '11371' : [16.5, 19.],
        '11372' : [15., 20.],
        '11373' : [15., 20.],
        '11374' : [15., 20.],
        '11375' : [15., 20.],
        '11376' : [15., 20.],
        '11377' : [15., 20.],
        '11378' : [15., 20.],
        '11379' : [15., 20.],
        '11380' : [15., 20.],
        '11381' : [15., 20.],
        '11382' : [15., 20.],
        '11383' : [15., 20.],
        '11384' : [15., 20.],
        '11385' : [15., 20.],
        '11386' : [15., 20.],
        '11387' : [15., 20.],
        '11388' : [15., 20.],
        '11389' : [15., 20.],
        '11390' : [15., 20.],
        '11391' : [15., 20.],
        '11392' : [15., 20.],
        '11393' : [15., 20.],
        '11394' : [15., 20.],
        '11395' : [15., 20.],
        '12932' : [15., 20.],
        '11397' : [15., 20.],
        '11398' : [15., 20.],
        '11399' : [15., 20.],
        '11400' : [15., 20.],
        '11401' : [15., 20.],
        '11402' : [15., 20.],
        '11403' : [15., 20.],
        '11404' : [15., 20.],
        '11405' : [15., 20.],
        '11406' : [15., 20.],
        '12943' : [15., 20.],
        '11408' : [15., 20.],
        '12945' : [15., 20.],
        '11412' : [15., 20.],
        '11413' : [15., 20.],
        '11419' : [15., 20.],
        '11420' : [15., 20.],
        '11421' : [15., 20.],
        '11422' : [15., 20.],
        '11423' : [15., 20.],
        '11424' : [15., 20.],
        '11425' : [15., 20.],
        '11426' : [15., 20.],
        '11427' : [15., 20.],
        '11428' : [15., 20.],
        '11429' : [15., 20.],
        '11430' : [15., 20.],
        '11431' : [15., 20.],
        '11432' : [15., 20.],
        '11433' : [15., 20.],
        '11434' : [15., 20.],
        '11435' : [15., 20.],
        '11436' : [15., 20.],
        '11437' : [15., 20.],
        '11438' : [15., 20.],
        '11440' : [15., 20.],
        '11441' : [15., 20.],
        '11442' : [15., 20.],
        '11443' : [15., 20.],
        '11444' : [15., 20.],
        '11446' : [15., 20.],
        '11447' : [15., 20.],
        '11448' : [15., 20.],
        '12923' : [15., 20.],
        '10969' : [15., 20.],
        '10970' : [15., 20.],
        '10971' : [15., 20.],
        '10972' : [15., 20.],
        '10973' : [15., 20.],
        '10974' : [15., 20.],
        '10975' : [15., 20.],
        '10976' : [15., 20.],
        '10977' : [15., 20.],
        '10978' : [15., 19.],
        '10979' : [15., 20.],
        '10980' : [15., 20.],
        '10983' : [15., 20.],
        '10984' : [15., 20.],
        '10985' : [15., 20.],
        '10986' : [15., 20.],
        '10987' : [15., 20.],
        '10988' : [15., 20.],
        '12925' : [15., 20.],
        '12926' : [15., 20.],
        '12927' : [15., 20.],
        '12928' : [15., 20.],
        '12929' : [15., 20.],
        '12930' : [15., 20.],
        '12933' : [15., 20.],
        '12934' : [15., 20.],
        '12935' : [15., 20.],
        '12937' : [15., 20.],
        '12938' : [15., 20.],
        '12940' : [15., 20.],
        '12941' : [15., 20.],
        '12942' : [15., 20.],
        '12944' : [15., 20.],
        '12924' : [15., 20.],
        '12922' : [15., 20.],
        '11220' : [15., 20.],
        '11221' : [15., 20.],
        '11222' : [15., 20.],
        '11223' : [15., 20.],
        '11224' : [15., 20.],
        '11225' : [15., 20.],
        '11226' : [15., 20.],
        '11227' : [15., 20.],
        '11228' : [15., 20.],
        '11229' : [15., 20.],
        '11230' : [15., 20.],
        '11231' : [15., 20.],
        '11232' : [15., 20.],
        '11233' : [15., 20.],
        '11234' : [15., 20.],
        '11235' : [15., 20.],
        '11236' : [15., 20.],
        '11237' : [15., 20.],
        '11238' : [15., 20.],
        '11239' : [15., 20.],
        '11240' : [15., 20.],
        '11242' : [15., 20.],
        '11243' : [15., 20.],
        '11244' : [15., 18.6],
        '11245' : [15., 20.],
        '11246' : [15., 20.],
        '11247' : [15., 20.],
        '11248' : [15., 20.],
        '11249' : [15., 20.],
        '11250' : [15., 20.],
        '11251' : [15., 20.],
        '11252' : [15., 20.],
        '11253' : [15., 20.],
        '11256' : [15., 20.],
        '11257' : [15., 20.],
        '11258' : [15., 20.],
        '11259' : [15., 20.],
        '11260' : [15., 20.],
        '11261' : [15., 20.],
        '11262' : [15., 20.],
        '11263' : [15., 20.]
    }.get(x, ' ')

    if out == ' ':
        raise ValueError("\n\nmagnitude limits for tile: "+x+" not found")
        
    return out






