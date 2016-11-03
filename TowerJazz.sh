#! /bin/bash

export PATHA="TowerJazz"

# Alignment
#export ALIG=$PATHA/cosmic_000768_000000.root
#export ALIG=$PATHA/cosmic_032873_000000.root
#export ALIG=$PATHA/cosmic_033078_000000.root
#export ALIG=$PATHA/cosmic_033040_000000.root
export ALIG=$PATHA/cosmic_033115_000000.root

# Setpoint specific
#export DUT=$PATHA/DUT
##export RCE=$PATHA/Cosmic
#export RCE=$PATHA/Cosmic_dataset2
#export RCE_MASKED=$PATHA/Cosmic_dataset2-mask
#export RCE_PROCC=$PATHA/Cosmic_dataset2-process
#export DUT_PROCC=$PATHA/DUT-process

#export DUT=$PATHA/run_804_tj_W3R15_50um_6V
#export DUT=$PATHA/dut_run804
#export DUT=$PATHA/dut_run804_tightenReset
#export DUT=$PATHA/dut_run804_first3k
#export DUT=$PATHA/dut_run804_right
#export DUT=$PATHA/run_953_tj_W3R13_20um_3V_1DRS_newStyle
#export DUT=$PATHA/run_804_tj_W3R15_50um_6V_new
#export DUT=$PATHA/run_868_tj_W3R13_50um_6V_1DRS_new
#export DUT=$PATHA/tmp_804
#export DUT=$PATHA/test_lowpass
#export DUT=$PATHA/run_804_tj_W3R15_50um_6V_fft4
#export DUT=$PATHA/run_32877_tj_W3R13_50um_6V_3DRS
#export DUT=$PATHA/run_32878_tj_W3R13_50um_6V_3DRS
#export DUT=$PATHA/test_30k
#export DUT=$PATHA/run_32878_tj_W3R13_50um_6V_3DRS_fft10
#export DUT=$PATHA/run_33072_tj_W3R13_M64_6V_3DRS_fft
#export DUT=$PATHA/run_33043_tj_W3R13_50um_6V_3DRS_fft
#export DUT=$PATHA/run_33118_tj_W3R15_M129_6V_3DRS_40k
export DUT=$PATHA/run_33118_tj_W3R15_M129_6V_3DRS_fft10
#export DUT=$PATHA/run_32878_tj_W3R13_50um_6V_3DRS_fft11
#export DUT=$PATHA/run_33072_tj_W3R13_M64_6V_3DRS_fft10
#export DUT=$PATHA/run_33090_tj_W3R13_M61_6V_3DRS_fft10
#export DUT=$PATHA/run_33072_tj_W3R13_M64_6V_3DRS_fft
#export DUT=$PATHA/run_32900_tj_W3R13_50um_6V_3DRS
#export DUT=$PATHA/run_868_tj_W3R13_50um_6V_1DRS_fft4
#export DUT=$PATHA/run804_new5
#export DUT=$PATHA/run868_new5
#export DUT=$PATHA/run_868_tj_W3R13_50um_6V_1DRS_fft4
#export DUT=$PATHA/shifted
#export DUT=$PATHA/inserted
#export DUT=$PATHA/dut_run925
#export RCE=$PATHA/cosmic_000804_000000
#export RCE=$PATHA/cosmic_000946_000000 # alignment
#export RCE=$PATHA/cosmic_000925_000000
#export RCE=$PATHA/cosmic_000868_000000
#export RCE=$PATHA/cosmic_000804_000000
#export RCE=$PATHA/cosmic_032877_000000
#export RCE=$PATHA/cosmic_032878_000000
#export RCE=$PATHA/cosmic_033072_000000
#export RCE=$PATHA/cosmic_033043_000000
export RCE=$PATHA/cosmic_033118_000000
#export RCE=$PATHA/cosmic_032878_000000
#export RCE=$PATHA/cosmic_033072_000000
#export RCE=$PATHA/cosmic_033080_000000
#export RCE=$PATHA/cosmic_033090_000000
#export RCE=$PATHA/cosmic_033072_000000
#export RCE=$PATHA/cosmic_032900_000000
export RCE_MASKED=${RCE}-mask
export RCE_PROCC=${RCE}-process
export DUT_PROCC=${DUT}-process

#export RCE=$PATHA/Cosmic
#export RCE=$PATHA/Cosmic
#export RCE_MASKED=$PATHA/Cosmic-mask
#export RCE_PROCC=$PATHA/Cosmic-process

#new
#export RCE=$PATHA/cosmic_000610
#export RCE_MASKED=$PATHA/cosmic_000610
#export RCE_PROCC=$PATHA/cosmic_000610-process

export DUTCONFIG=dutTowerJazz.cfg
#export DUTCONFIG=dutTowerJazz20um.cfg
#export DUTCONFIG=dutTowerJazz25um.cfg
#export DUTCONFIG=dutTowerJazz30um.cfg
#export DUTCONFIG=dutTowerJazz.cfg

#export TELECONF=configs/reforig_TowerJazz.cfg
export TELECONF=configs/reforig_TowerJazz5b.cfg
#export TELECONF=configs/reforig_TowerJazz5c.cfg

echo "TowerJazz Efficiency Analysis --- 0 TILT"

echo "-------------------Apply mask"
#./Judith -c applyMask -i $RCE.root -o ${RCE_MASKED}.root  -r $TELECONF -t configs/globalorig_TowerJazz.cfg #-n 50000

echo "-------------------syncronization"
#./Judith -c synchronizeRMS -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -s 0 #-n 80000
./Judith -c synchronizeBunchXing -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -s 0 #-n 40000
#./Judith -c synchronizeBunchXing -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -s 0 -n 30000
#./Judith -c synchronize -i ${RCE_MASKED}.root -o ${RCE_MASKED}_sync.root -I $DUT.root -O ${DUT}_sync.root -r $TELECONF -d configs/$DUTCONFIG -t  configs/globalorig_TowerJazz.cfg -s 0 # -n 4000

echo "-------------------CoarseAlign telescope"
./Judith -c coarseAlign -i $ALIG -r $TELECONF -t configs/globalorig_TowerJazz.cfg  #-n 100000
echo "-------------------FineAlign telescope"
#./Judith -c fineAlign -i $ALIG -r $TELECONF -t configs/globalorig_TowerJazz.cfg -R ${ALIG}-alignment-result.root # -n 50000
#
## change to the sync names
export RCE_MASKED=${RCE_MASKED}_sync
export DUT=${DUT}_sync
echo "-------------------CoarseAlign dut"
./Judith -c coarseAlignDUT -i ${RCE_MASKED}.root -I ${DUT}.root -r $TELECONF -t configs/globalorig.cfg -d configs/$DUTCONFIG  #-n 10000 #--eventOffset 50000
#echo "-------------------FineAlign dut"
./Judith -c fineAlignDUT -i ${RCE_MASKED}.root -I ${DUT}.root -t configs/globalorig.cfg -r $TELECONF -d configs/$DUTCONFIG -R DUT-alignment-result.root #-n 50000

#export RCE_MASKED=$PATHA/cosmic_033040_000000-mask
#export RCE_PROCC=$PATHA/cosmic_033040_000000-process
#export RCE=$PATHA/cosmic_033040_000000-mask
#echo "-------------------process telescope"
./Judith -c process -i ${RCE_MASKED}.root -o ${RCE_PROCC}.root -r $TELECONF -t configs/globalorig_TowerJazz.cfg -R ${RCE}-proccess-result.root #-n 50000
#echo "-------------------process dut"
./Judith -c process -i ${DUT}.root -o ${DUT_PROCC}.root -r configs/$DUTCONFIG -t configs/globalorig.cfg -R ${DUT}-procress-result.root 
##-n 50000
#echo "-------------------Analysis telescope"
./Judith -c analysis -i ${RCE_PROCC}.root -r $TELECONF -t configs/globalorig.cfg -R ${RCE}-analysis-result.root #-n 50000
#echo "-------------------Analysis dut"
./Judith -c analysisDUT -i ${RCE_PROCC}.root -I ${DUT_PROCC}.root -r $TELECONF -d configs/$DUTCONFIG -t configs/globalorig_TowerJazz.cfg -R ${DUT}-analysis-result.root #-n 28000 #-s 40000 #-n 100000 #--eventOffset 50000  #-n 60000

