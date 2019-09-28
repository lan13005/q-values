sed -i "s@varNameSet@cosTheta_X_cms;phi_X_cms;cosTheta_eta_gjs;phi_eta_gjs;cosThetaHighestEphotonIneta_gjs;cosThetaHighestEphotonInpi0_cms@g" run.sh
./run.sh
sed -i "s@cosTheta_X_cms;phi_X_cms;cosTheta_eta_gjs;phi_eta_gjs;cosThetaHighestEphotonIneta_gjs;cosThetaHighestEphotonInpi0_cms@varNameSet@g" run.sh
