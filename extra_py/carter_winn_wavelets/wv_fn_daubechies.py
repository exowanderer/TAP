# Id: //depot/idl/IDL_71/idldir/lib/wavelet/source/wv_fn_daubechies.pro  # 1 \
#
# Copyright (c) 1999-2009, ITT Visual Information Solutions. All
#       rights reserved. Unauthorized reproduction is prohibited.
#+
# NAME:
#    WV_FN_DAUBECHIES
#
# PURPOSE:
#    This def  returns the Daubechies wavelet coefficients.
#
# CALLING SEQUENCE:
#    info = WV_FN_DAUBECHIES( Order, Scaling, Wavelet, Ioff, Joff)
#
# INPUTS:
#    Order: This is the order number for the Daubechies wavelet.
#           Order=2 is the "db2" wavelet, and has 2 vanishing moments
#           and 4 coefficients.
#
# OUTPUTS:
#    Scaling: A vector of the scaling (father) coefficients
#
#    Wavelet: A vector of the wavelet (mother) coefficients
#
#    Ioff: The offset index used to center the Scaling support
#
#    Joff: The offset index used to center the Wavelet support
#
# KEYWORD PARAMETERS:
#    None.
#
# RETURN VALUE:
#    Returns a structure with the following information:
#          (this is an example for order=2)
#       info = {family:'Daubechies', \    # name of wavelet family
#               order_name:'Order', \     # term used for "order"
#               order_range:[2,24,2], \   # valid range [first,last,default]
#               order:2, \                # order number
#               discrete:1, \             # 0=continuous, 1=discrete
#               orthogonal:1, \           # 0=nonorthogonal, 1=orthogonal
#               symmetric:0, \            # 0=asymmetric, 1=symmetric
#               support:3, \              # support width
#               moments:2, \              # # of vanishing moments
#               regularity:0.550}         # # of continuous derivatives
#
# REFERENCE:
#    Daubechies, I., 1992: Ten Lectures on Wavelets, SIAM, p. 195.
#       Daubechies has orders 2-10, although note that Daubechies
#       has multiplied by Sqrt(2).
#
#    Orders 11-24 are from &lt#http://www.isds.duke.edu/~brani/filters.html&gt#
#
# MODifICATION HISTORY:
#    Written by: Chris Torrence, 1999
#-


def info_message(message):
    print(f'[INFO] {message}')

#----------------------------------------------------------------
# Daubechies orthogonal asymmetric wavelet coefficients


def wv_fn_daubechies_coeff(order):
    # COMPILE_OPT strictarr, hidden
    if order == 1:
        coeff = [1, 1] / np.sqrt(2)  # same as Haar
    elif order == 2:
        coeff = [1 + np.sqrt(3), 3 + np.sqrt(3),
                 3 - np.sqrt(3), 1 - np.sqrt(3)] / (4 * np.sqrt(2))
    elif order == 3:
        sq10 = np.sqrt(10)
        sq5210 = SQRT(5 + 2 * np.sqrt(10))
        coeff = [1 + sq10 + sq5210, 5 + sq10 + 3 * sq5210,
                 10 - 2 * sq10 + 2 * sq5210, 10 - 2 * sq10 - 2 * sq5210,
                 5 + sq10 - 3 * sq5210, 1 + sq10 - sq5210] / (16 * np.sqrt(2))
    elif order == 4:
        coeff = [
            0.2303778133088965,
            0.7148465705529158,
            0.630880767929859,
            -0.02798376941686011,
            -0.1870348117190932,
            0.0308413818355608,
            0.03288301166688522,
            -0.01059740178506904]
    elif order == 5:
        coeff = [
            0.1601023979741924,
            0.6038292697971881,
            0.7243085284377715,
            0.1384281459013217,
            -0.2422948870663802,
            -0.03224486958463778,
            0.07757149384004565,
            -0.006241490212798174,
            -0.01258075199908194,
            0.003335725285473757]
    elif order == 6:
        coeff = [
            0.11154074335011,
            0.4946238903984554,
            0.7511339080210982,
            0.315250351709197,
            -0.2262646939654429,
            -0.1297668675672638,
            0.0975016055873231,
            0.02752286553030565,
            -0.03158203931748625,
            0.0005538422011615105,
            0.004777257510945544,
            -0.001077301085308486]
    elif order == 7:
        coeff = [
            0.07785205408500813,
            0.3965393194819123,
            0.7291320908462274,
            0.4697822874051917,
            -0.1439060039285563,
            -0.2240361849938672,
            0.07130921926683042,
            0.080612609151082,
            -0.03802993693501439,
            -0.016574541630667,
            0.01255099855609955,
            0.0004295779729213739,
            -0.001801640704047446,
            0.0003537137999745171]
    elif order == 8:
        coeff = [
            0.05441584224310704,
            0.3128715909143165,
            0.6756307362973218,
            0.5853546836542239,
            -0.01582910525637238,
            -0.2840155429615815,
            0.0004724845739030209,
            0.1287474266204823,
            -0.01736930100181088,
            -0.04408825393079791,
            0.01398102791739956,
            0.00874609404740648,
            -0.004870352993451852,
            -0.000391740373376942,
            0.0006754494064506183,
            -0.0001174767841247786]
    elif order == 9:
        coeff = [
            0.03807794736388813,
            0.2438346746126514,
            0.6048231236902548,
            0.6572880780514298,
            0.1331973858249681,
            -0.2932737832793372,
            -0.0968407832230689,
            0.148540749338104,
            0.03072568147931585,
            -0.06763282906135907,
            0.0002509471148277948,
            0.02236166212368439,
            -0.004723204757752752,
            -0.004281503682464633,
            0.001847646883056686,
            0.0002303857635232296,
            -0.0002519631889427889,
            0.00003934732031628112]
    elif order == 10:
        coeff = [
            0.02667005790054869,
            0.188176800077648,
            0.527201188931628,
            0.6884590394535462,
            0.2811723436606982,
            -0.2498464243271048,
            -0.1959462743773243,
            0.127369340335694,
            0.0930573646035142,
            -0.07139414716638016,
            -0.0294575368218849,
            0.03321267405931551,
            0.003606553566951515,
            -0.0107331754833277,
            0.001395351747051327,
            0.001992405295184184,
            -0.0006858566949593225,
            -0.0001164668551292262,
            0.0000935886703200315,
            -0.00001326420289451403]
    elif order == 11:
        coeff = [
            0.01869429776144806,
            0.1440670211504498,
            0.4498997643555165,
            0.6856867749154562,
            0.4119643689476272,
            -0.1622752450269621,
            -0.2742308468172826,
            0.06604358819685894,
            0.1498120124663909,
            -0.04647995511648684,
            -0.06643878569486228,
            0.03133509021904213,
            0.02084090436017028,
            -0.01536482090617611,
            -0.003340858873009247,
            0.0049284176560525,
            -0.0003085928588149355,
            -0.00089302325066525,
            0.0002491525235524301,
            0.00005443907469928305,
            -0.00003463498418694142,
            0.000004494274277230458]
    elif order == 12:
        coeff = [
            0.01311225795736534,
            0.1095662728222715,
            0.3773551352176745,
            0.657198722584349,
            0.5158864784293156,
            -0.04476388565908393,
            -0.3161784537592869,
            -0.02377925725693821,
            0.1824786059298069,
            0.00535956967427179,
            -0.0964321200976865,
            0.0108491302560784,
            0.04154627749559747,
            -0.01221864906995923,
            -0.01284082519846823,
            0.00671149900888981,
            0.002248607241020708,
            -0.002179503618657147,
            0.000006545128213682533,
            0.0003886530628261407,
            -0.0000885041092094801,
            -0.00002424154575734139,
            0.00001277695221955214,
            -0.000001529071758089919]
    elif order == 13:
        coeff = [
            0.00920213353936357,
            0.082861243876398,
            0.3119963221728867,
            0.6110558511805082,
            0.5888895704451372,
            0.0869857261666496,
            -0.314972907739053,
            -0.124576730762086,
            0.1794760794355785,
            0.07294893365742099,
            -0.1058076181950538,
            -0.02648840647689916,
            0.05613947710301562,
            0.002379972253836755,
            -0.02383142071161908,
            0.003923941449079961,
            0.007255589402002825,
            -0.002761911234808676,
            -0.001315673911943637,
            0.000932326130928484,
            0.00004925152513188404,
            -0.0001651289885636495,
            0.00003067853758174376,
            0.00001044193057207212,
            -0.000004700416479607929,
            0.0000005220035098765021]
    elif order == 14:
        coeff = [
            0.006461153459818989,
            0.0623647588469322,
            0.2548502677833766,
            0.5543056179241174,
            0.6311878490950694,
            0.2186706877760189,
            -0.2716885522429336,
            -0.2180335299738394,
            0.138395213856541,
            0.1399890165735457,
            -0.0867484115685856,
            -0.07154895550625034,
            0.05523712625188016,
            0.02698140830446938,
            -0.0301853515397028,
            -0.005615049530747707,
            0.01278949326524909,
            -0.000746218989436958,
            -0.003849638867994312,
            0.001061691085418039,
            0.0007080211541344865,
            -0.0003868319473184179,
            -0.00004177724577935138,
            0.00006875504251988474,
            -0.00001033720918460207,
            -0.000004389704901652653,
            0.000001724994675254821,
            -0.000000178713996820958]
    elif order == 15:
        coeff = [
            0.004538537356680069,
            0.0467433948433292,
            0.2060238637760462,
            0.4926317712332494,
            0.6458131398235114,
            0.339002535383428,
            -0.19320413905893,
            -0.2888825960016258,
            0.06528295291444258,
            0.1901467139017971,
            -0.03966617641454303,
            -0.1111209358626346,
            0.03387714389352461,
            0.05478055052762776,
            -0.0257670072911817,
            -0.02081005014572826,
            0.01508391800773139,
            0.005101000354434229,
            -0.006487734552531616,
            -0.0002417564910950625,
            0.001943323977748212,
            -0.0003734823537271217,
            -0.0003595652439869339,
            0.0001558964896924794,
            0.00002579269911910246,
            -0.00002813329623232866,
            0.000003362987176654478,
            0.000001811270405641324,
            -0.0000006316882317817563,
            0.00000006133359905269254]
    elif order == 16:
        coeff = [
            0.003189220905181802,
            0.0349077141074775,
            0.1650642824989111,
            0.4303127204089899,
            0.6373563289234388,
            0.4402902557886062,
            -0.0897510867287953,
            -0.3270633068118058,
            -0.02791820715372535,
            0.2111906930487478,
            0.02734026408611786,
            -0.1323883043443139,
            -0.00623972263724492,
            0.07592423555847598,
            -0.00758897425298305,
            -0.03688839741760147,
            0.01029765955546528,
            0.01399376876290007,
            -0.006990014507518413,
            -0.003644279596729619,
            0.003128023357662664,
            0.000407896978913364,
            -0.000941021742187743,
            0.000114241519113091,
            0.0001747872440135933,
            -0.00006103596571228747,
            -0.00001394566888488284,
            0.00001133660857799308,
            -0.000001043571333041443,
            -0.0000007363656730469882,
            0.0000002308784069376313,
            -0.00000002109339613774453]
    elif order == 17:
        coeff = [
            0.002241806968367765,
            0.02598539333038641,
            0.1312149014643511,
            0.3703507191428474,
            0.6109966080619875,
            0.5183157592365552,
            0.02731497388861195,
            -0.3283207398752789,
            -0.1265997478695799,
            0.1973105883690036,
            0.1011354893285621,
            -0.1268156885448092,
            -0.05709141812622551,
            0.081105985705437,
            0.02231233608959475,
            -0.04692243752178137,
            -0.003270955473782776,
            0.02273367623263168,
            -0.003042989911563062,
            -0.00860292137975392,
            0.002967996640915282,
            0.002301205207197428,
            -0.001436845280352317,
            -0.0003281325149411173,
            0.0004394654201169656,
            -0.00002561010931458864,
            -0.0000820480308801988,
            0.00002318681330990614,
            0.000006990600842366534,
            -0.000004505942411707292,
            0.0000003016549532645506,
            0.0000002957700881589635,
            -0.0000000842394830828037,
            0.000000007267492843919008]
    elif order == 18:
        coeff = [
            0.001576310332632241,
            0.01928853309434481,
            0.1035884729715391,
            0.3146789620466176,
            0.571826841995251,
            0.5718016803655575,
            0.147223099399332,
            -0.2936540837163994,
            -0.2164809618743174,
            0.1495339814252923,
            0.1670813196471977,
            -0.0923318969776604,
            -0.1067522571200224,
            0.0648872212223416,
            0.05705125157931265,
            -0.04452614611490133,
            -0.02373321210978654,
            0.02667070832113655,
            0.006262168357742094,
            -0.01305148206344844,
            0.0001186301071328846,
            0.004943344018360076,
            -0.001118732786346494,
            -0.001340596411265555,
            0.0006284657384942994,
            0.0002135815764103265,
            -0.000198648570821057,
            -0.000000153591634265962,
            0.00003741238184339052,
            -0.00000852060341054129,
            -0.00000333263477007513,
            0.00000176871313748643,
            -0.00000007691633640217469,
            -0.0000001176098869880653,
            0.00000003068836137122469,
            -0.000000002507934683892356]
    elif order == 19:
        coeff = [
            0.001108669779715294,
            0.01428109865333334,
            0.081278114333354,
            0.2643884347822977,
            0.5244363819574067,
            0.6017045501513535,
            0.2608949440110274,
            -0.2280914100170829,
            -0.285838641929714,
            0.07465227262054114,
            0.2123497512548378,
            -0.03351853842979753,
            -0.1427856935054576,
            0.02758435493215239,
            0.0869067594236619,
            -0.02650123589611068,
            -0.04567422669495623,
            0.02162376812192859,
            0.01937555029280247,
            -0.01398838901012597,
            -0.00586692239134182,
            0.007040747519198927,
            0.0007689543646753964,
            -0.002687551858597481,
            0.0003418086639330359,
            0.0007358025360798398,
            -0.0002606761416764582,
            -0.0001246007941078683,
            0.0000871127066319985,
            0.000005105950548947162,
            -0.00001664017665533139,
            0.000003010964385934741,
            0.000001531931507655374,
            -0.0000006862755810090276,
            0.00000001447088339408005,
            0.00000004636937873589416,
            -0.000000011164020912898,
            0.000000000866684902796269]
    elif order == 20:
        coeff = [
            0.0007799530020084384,
            0.0105493864101072,
            0.06342373157542249,
            0.2199419467839922,
            0.4726958375631425,
            0.6104928215175741,
            0.3615021297395791,
            -0.139211825416023,
            -0.3267863905078842,
            -0.01672694530514085,
            0.2282909876975237,
            0.03985032729018178,
            -0.1554585361790331,
            -0.02471674917392653,
            0.1022916746204368,
            0.005632268726873665,
            -0.06172283526148656,
            0.005874682288534986,
            0.03229427583633914,
            -0.00878931595226129,
            -0.01381051445886118,
            0.006721621652169426,
            0.004420538864131319,
            -0.003581491222634283,
            -0.00083156152944895,
            0.001392558453825609,
            -0.00005349753868856166,
            -0.0003851044297986765,
            0.0001015328014373285,
            0.00006774275277093538,
            -0.00003710583043522718,
            -0.000004376140493506968,
            0.000007241242222701708,
            -0.000001011993125412585,
            -0.0000006847073928591012,
            0.0000002633921999175421,
            0.0000000002014328820034285,
            -0.0000000181484172957345,
            0.000000004056123630675098,
            -0.0000000002998833944499773]
    elif order == 21:
        coeff = [
            0.0005488240399453808,
            0.007776660464348811,
            0.04924790475876491,
            0.1813601028599902,
            0.419688998145241,
            0.6015074510688103,
            0.4445910837993439,
            -0.03572381948901234,
            -0.33566575122537,
            -0.1123978514710653,
            0.2115648260162405,
            0.1152333439473735,
            -0.1399410472763452,
            -0.08177625782428998,
            0.09660066710664022,
            0.04572352417673011,
            -0.06497770623152748,
            -0.01865389796875268,
            0.03972696757220106,
            0.003357765554657301,
            -0.02089211624987374,
            0.002403482102825579,
            0.008988852342563074,
            -0.002891344156898007,
            -0.002958382842307337,
            0.001716612683276365,
            0.0006394203289590759,
            -0.0006906733219030776,
            -0.00003196410553726866,
            0.0001936652571660039,
            -0.0000363553295677002,
            -0.00003499676704742804,
            0.00001535487521020741,
            2.79033850314008 - 6,
            -3.090027001911602 - 6,
            3.16610662424439 - 7,
            2.99214595113828 - 7,
            -1.000404119487493 - 7,
            -2.254019869522092 - 9,
            7.058055911572644 - 9,
            -1.471958939283684 - 9,
            1.038808947669207 - 10]
    elif order == 22:
        coeff = [
            0.0003862673246197253,
            0.005721914066631049,
            0.03807032187691932,
            0.1483689789282081,
            0.3677320057234413,
            0.5784372354311235,
            0.5079033273631367,
            0.07372115020105462,
            -0.3127333476121842,
            -0.2005720141344328,
            0.1640948426591233,
            0.179974931810035,
            -0.0971123372197599,
            -0.1317696149504392,
            0.06807740848784511,
            0.08455839833964807,
            -0.05136497255398131,
            -0.04653131832736136,
            0.03697137276735332,
            0.02058693268949487,
            -0.02348031395539096,
            -0.006213835918293782,
            0.01256489065516637,
            0.0003001305020824184,
            -0.005455761185358356,
            0.001044278408986017,
            0.001827032986409597,
            -0.000770702101944467,
            -0.0004237923063271874,
            0.0003286138886837352,
            0.0000434593692542139,
            -0.00009405347080647135,
            0.00001137454223403893,
            0.00001737397675279249,
            -6.166816318076451 - 6,
            -1.565197277819435 - 6,
            1.295199441207159 - 6,
            -8.78003044824892 - 8,
            -1.283352833826032 - 7,
            3.761280659022215 - 8,
            1.680187679020641 - 9,
            -2.729659356918619 - 9,
            5.33601149622179 - 10,
            -3.60216327759258 - 11]
    elif order == 23:
        coeff = [
            0.0002719278182602901,
            0.004203109552950134,
            0.02931247643736339,
            0.1205254471036576,
            0.3184759568589838,
            0.5449708209347766,
            0.5510501337055957,
            0.1813841378320262,
            -0.2614398761995617,
            -0.2714429864972958,
            0.0921245749243952,
            0.2235864349031235,
            -0.03304774793732929,
            -0.164030308293076,
            0.02028436820991752,
            0.1123069840244809,
            -0.0211292480280753,
            -0.07021415427385447,
            0.02176834655240395,
            0.03849895908078205,
            -0.01852549112315692,
            -0.01753870400805271,
            0.01275326613768589,
            0.006032371159860696,
            -0.00707603267773538,
            -0.001134947880346942,
            0.003123184807392083,
            -0.000246537026777104,
            -0.001061334361043996,
            0.000319454992361999,
            0.0002567865998070605,
            -0.0001500373078747796,
            -0.00003379188332733358,
            0.00004426515179248939,
            -2.635561787093299 - 6,
            -8.348692795439998 - 6,
            2.397822036092728 - 6,
            8.148343038224153 - 7,
            -5.339546450998099 - 7,
            1.853340229309923 - 8,
            5.418084825798256 - 8,
            -1.400079829615052 - 8,
            -9.473736128438874 - 10,
            1.050551729533758 - 9,
            -1.93260193304542 - 10,
            1.250331739337031 - 11]
    elif order == 24:
        coeff = [
            0.0001914240079776934,
            0.003081894336144903,
            0.02248099723913652,
            0.09725657409395711,
            0.272893661713225,
            0.5043448957614517,
            0.5749146829767083,
            0.2809851510053765,
            -0.1872418464658568,
            -0.3179111538203686,
            0.004781510762825361,
            0.2392258659829295,
            0.042531243536347,
            -0.1711600617797226,
            -0.03877318682438014,
            0.1210092088290207,
            0.02098022912439134,
            -0.08215538086453539,
            -0.004578395730450242,
            0.05129798128535279,
            -0.004944235600686442,
            -0.02821125709939177,
            0.007661004281903052,
            0.01304905186620713,
            -0.006290964935213451,
            -0.004746267936383896,
            0.00373576397589871,
            0.001153694353296646,
            -0.001696334910033699,
            -0.00004416435334971148,
            0.0005860851561798487,
            -0.000118113728929818,
            -0.0001459980983446589,
            0.00006558881863639525,
            0.00002183096348720674,
            -0.00002022741617379432,
            1.337052417608915 - 8,
            3.900821594914755 - 6,
            -8.979550384702172 - 7,
            -4.032228084773544 - 7,
            2.166180932866001 - 7,
            -5.054643465620961 - 10,
            -2.255577015054618 - 8,
            5.157391468496204 - 9,
            4.748066278754132 - 10,
            -4.024365393060184 - 10,
            6.991284124010881 - 11,
            -4.342457865150871 - 12]
    else:
        coeff = -1

    return coeff


#----------------------------------------------------------------
def wv_fn_daubechies(order, scaling, wavelet, ioff, joff):
    # COMPILE_OPT strictarr
    # ON_ERROR, 2  # return to caller

    # defaults
    order_range = [1, 24, 2]  # [first,last,default]
    if(len(order) < 1):
        order = order_range[2]  # default

    # check for invalid Order
    # order = FIX(order + 1e-5)
    order = int(order + 1e-5)
    if order < order_range[0] or order > order_range[1]:
        info_message('Order out of range, reverting to default...')
        order = order_range[2]  # default

    # Regularity = # of continuous derivatives
    # orders  0-4 from Daubechies (1992) p. 239
    # orders 5-24 from Taswell (1998) p. 36 "DROMA"
    regularity24 = [0, 0, 0.550, 1.088,   # order 0-3
                    1.618, 1.969, 2.189, 2.460,    # 4-7
                    2.761, 3.074, 3.381, 3.603,    # 8-11
                    3.833, 4.073, 4.317, 4.558,    # 12-15
                    4.791, 5.014, 5.239, 5.465,    # 16-19
                    5.691, 5.916, 6.138, 6.360,    # 20-23
                    6.581]  # 24
    if order <= 24:
        regularity = regularity24[order]

    else:  # orders > 24 from Daubechies (1992) p. 226
        regularity = (0.2 * order)
        regularity = regularity > regularity24[24]  # CT assumption

    # construct Info structure
    info = {'family': 'Daubechies',
            'order_name': 'Order',
            'order_range': order_range,
            'order': order,
            'discrete': 1,
            'orthogonal': 1,
            'symmetric': 0,
            'support': int(order * 2 - 1),
            'moments': int(order),
            'regularity': regularity}

    # if(N_PARAMS() < 1) THEN RETURN, info  # don't bother with rest
    if order is None:
        return info

    # choose scaling coefficients
    scaling = wv_fn_daubechies(order, scaling, wavelet, ioff, joff)
    # scaling = wv_fn_daubechies_coeff(order)

    # construct wavelet coefficients & offsets
    n = len(scaling)
    wavelet = scaling[::-1] * ((-1)**range(n))
    ioff = int(-n // 2 + 2)  # offset for scaling
    joff = ioff      # offset for wavelet

    return info, n, wavelet, ioff, joff