#/bin/bash 

Rscript Code/Settings/1_getData.R

if [ "$1" == "1" ]; then 
	Rscript Code/Options/2_1_SubsettingTumourSubtypes.R $2 $3

elif [ "$1" == "2" ]; then
	Rscript Code/Options/2_2_SubsettingNodeNumber.R $2 $3

elif [ "$1" == "3" ]; then
	Rscript Code/Options/2_3_subsettingAges.R $2 $3

elif [ "$1" == "4" ]; then
	Rscript Code/Options/2_4_subsettingGrade.R

elif [ "$1" == "5" ]; then
	Rscript Code/Options/2_5_SubsettingNodes.R

elif [ "$1" == "6" ]; then
	Rscript Code/Options/2_6_SubsettingERStatus.R

else
	Rscript Code/Options/2_7_SubsettingPRStatus.R
fi

Rscript Code/Options/3_DEA.R