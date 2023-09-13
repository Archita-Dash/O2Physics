# fulljets
This repo is for my correlation studies of dileptons inside of fulljets

# AO2D.root
This is Run3_pass4 data file downloaded from https://alimonitor.cern.ch/catalogue/index.jsp?path=%2Falice%2Fdata%2F2022%2FLHC22m%2F523308%2Fapass4_tpc_v1%2F0150#/alice/data/2022/LHC22m/523308/apass4_tpc_v1/0150/o2_ctf_run00523308_orbit0363126684_tf0000344007_epn038/001 . It's a pass 4 data and hence is good!
Steps to download the Run3 data files could be found here: https://github.com/zchochul/AliceO2

##As mentioned below problems with the above AO2D file. On 13th Sept. 2023, I downloaded the following latest AO2D file from this year which is pass1
/alice/data/2023/LHC23zk/539218/apass1/0210/o2_ctf_run00539218_orbit0056460544_tf0000000048_epn213/001

# How to run
1. Enter O2Physics/latest environment
2. bash <runscript.sh> <inputdata-file> <dpl-config.json>

#How to run the full jet trigger task
1. Enter O2Physics/latest environment
2. bash runfulljet_script.sh AO2D.root dpl-config.json		// on my office desktop: cd /data/alice/fulljets
3. The above cmd produces the AnalysisResults.root file		//Tried with the above AO2D.root file from LHC22m but somehow the "full-jet-filter" task is empty
						// Markus suggested to try with a data set later than LHC22o since there were problems with "m" and "o" datasets
  
