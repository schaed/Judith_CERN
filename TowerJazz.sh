#! /bin/bash

export PATH="TowerJazz"

# Alignment
#export ALIG=$PATH/cosmic_000768_000000.root
export ALIG=$PATH/cosmic_032873_000000.root

# Setpoint specific
#export DUT=$PATH/DUT
##export RCE=$PATH/Cosmic
#export RCE=$PATH/Cosmic_dataset2
#export RCE_MASKED=$PATH/Cosmic_dataset2-mask
#export RCE_PROCC=$PATH/Cosmic_dataset2-process
#export DUT_PROCC=$PATH/DUT-process

#export DUT=$PATH/run_804_tj_W3R15_50um_6V
#export DUT=$PATH/dut_run804
#export DUT=$PATH/dut_run804_tightenReset
#export DUT=$PATH/dut_run804_first3k
#export DUT=$PATH/dut_run804_right
#export DUT=$PATH/run_953_tj_W3R13_20um_3V_1DRS_newStyle
#export DUT=$PATH/run_804_tj_W3R15_50um_6V_new
#export DUT=$PATH/run_868_tj_W3R13_50um_6V_1DRS_new
#export DUT=$PATH/tmp_804
#export DUT=$PATH/test_lowpass
#export DUT=$PATH/run_804_tj_W3R15_50um_6V_fft4
#export DUT=$PATH/run_32877_tj_W3R13_50um_6V_3DRS
#export DUT=$PATH/run_32878_tj_W3R13_50um_6V_3DRS
#export DUT=$PATH/test_30k
export DUT=$PATH/test_30k_3
#export DUT=$PATH/run_32900_tj_W3R13_50um_6V_3DRS
#export DUT=$PATH/run_868_tj_W3R13_50um_6V_1DRS_fft4
#export DUT=$PATH/run804_new5
#export DUT=$PATH/run868_new5
#export DUT=$PATH/run_868_tj_W3R13_50um_6V_1DRS_fft4
#export DUT=$PATH/shifted
#export DUT=$PATH/inserted
#export DUT=$PATH/dut_run925
#export RCE=$PATH/cosmic_000804_000000
#export RCE=$PATH/cosmic_000946_000000 # alignment
#export RCE=$PATH/cosmic_000925_000000
#export RCE=$PATH/cosmic_000868_000000
#export RCE=$PATH/cosmic_000804_000000
#export RCE=$PATH/cosmic_032877_000000
export RCE=$PATH/cosmic_032878_000000
#export RCE=$PATH/cosmic_032900_000000
export RCE_MASKED=${RCE}-mask
export RCE_PROCC=${RCE}-process
export DUT_PROCC=${DUT}-process

#export RCE=$PATH/Cosmic
#export RCE=$PATH/Cosmic
#export RCE_MASKED=$PATH/Cosmic-mask
#export RCE_PROCC=$PATH/Cosmic-process

#new
#export RCE=$PATH/cosmic_000610
#export RCE_MASKED=$PATH/cosmic_000610
#export RCE_PROCC=$PATH/cosmic_000610-process

#export DUTCONFIG=dutTowerJazz.cfg
#export DUTCONFIG=dutTowerJazz20um.cfg
export DUTCONFIG=dutTowerJazz.cfg

#export TELECONF=configs/reforig_TowerJazz.cfg
export TELECONF=configs/reforig_TowerJazz5b.cfg

echo "TowerJazz Efficiency Analysis --- 0 TILT"

echo "-------------------Apply mask"
#./Judith -c applyMask -i $RCE.root -o ${RCE_MASKED}.root  -r $TELECONF -t configs/globalorig_TowerJazz.cfg #-n 50000

echo "-------------------syncronization"
#./Judith -c synchronizeRMS -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -s 0 #-n 80000
#./Judith -c synchronizeBunchXing -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -s 0 #-n 10000
./Judith -c synchronizeBunchXing -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -s 0 -n 30000
#./Judith -c synchronize -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t  configs/globalorig_TowerJazz.cfg -s 0 # -n 4000

echo "-------------------CoarseAlign telescope"
#./Judith -c coarseAlign -i $ALIG -r $TELECONF -t configs/globalorig_TowerJazz.cfg  #-n 100000
echo "-------------------FineAlign telescope"
#./Judith -c fineAlign -i $ALIG -r $TELECONF -t configs/globalorig_TowerJazz.cfg -R ${ALIG}-alignment-result.root # -n 50000
#
## change to the sync names
export RCE_MASKED=${RCE_MASKED}_sync
export DUT=${DUT}_sync
echo "-------------------CoarseAlign dut"
#./Judith -c coarseAlignDUT -i ${RCE_MASKED}.root -I ${DUT}.root -r $TELECONF -t configs/globalorig.cfg -d configs/$DUTCONFIG
#echo "-------------------FineAlign dut"
#./Judith -c fineAlignDUT -i ${RCE_MASKED}.root -I ${DUT}.root -t configs/globalorig.cfg -r $TELECONF -d configs/$DUTCONFIG -R DUT-alignment-result.root -n 50000


echo "-------------------process telescope"
./Judith -c process -i ${RCE_MASKED}.root -o ${RCE_PROCC}.root -r $TELECONF -t configs/globalorig_TowerJazz.cfg -R ${RCE}-proccess-result.root #-n 50000
echo "-------------------process dut"
./Judith -c process -i ${DUT}.root -o ${DUT_PROCC}.root -r configs/$DUTCONFIG -t configs/globalorig.cfg -R ${DUT}-procress-result.root
#-n 50000
echo "-------------------Analysis telescope"
./Judith -c analysis -i ${RCE_PROCC}.root -r $TELECONF -t configs/globalorig.cfg -R ${RCE}-analysis-result.root #-n 50000
echo "-------------------Analysis dut"
./Judith -c analysisDUT -i ${RCE_PROCC}.root -I ${DUT_PROCC}.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -R ${DUT}-analysis-result.root #-n 60000

