
//for compact_Lhe_dip_sol_ell.xml
//-------------------------------------------------------------------------------------------------
{gROOT->Reset();
   
/***********************************calculations of silicon area**************************************/
//-----------------------hcalplug_fwd--------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

float VxBar_Env= 0.5;
float CentralBeamPipe_rmax=2.4;
float Distance_VXDLayer= 6;
float VertexBarrel_r0=CentralBeamPipe_rmax+ VxBar_Env + 0.6;
float VertexBarrel_r1=VertexBarrel_r0 + VxBar_Env + Distance_VXDLayer;
float VertexBarrel_r2=VertexBarrel_r1 + VxBar_Env + Distance_VXDLayer;
 
float Distance_SITLayer=6.0;

int CaloSides=12;
float VXD_Ell_Max_r=VertexBarrel_r0;
float Radius_SITLayer0=VXD_Ell_Max_r + 5.5;
float Radius_SITLayer1=Radius_SITLayer0 + Distance_SITLayer;
float Radius_SITLayer2=Radius_SITLayer1 + Distance_SITLayer;
float Radius_SITLayer3=Radius_SITLayer2 + Distance_SITLayer;
float Radius_SITLayer4=Radius_SITLayer3 + Distance_SITLayer;

float Diff_Radius_SOTLayer= 6.;
float Distance_SOTLayer = 15.;
float Radius_SOTLayer0 = Radius_SITLayer2 + Diff_Radius_SOTLayer;
float Radius_SOTLayer1 = Radius_SOTLayer0 + Distance_SOTLayer;
float Radius_SOTLayer2 = Radius_SOTLayer1 + Distance_SOTLayer;

float EcalBarrel_rmin=((Radius_SOTLayer2 + 6.0) / (TMath::Cos(3.14159/CaloSides)) );
//float angle= TMath::Cos(3.14159/CaloSides);
//cout << " angle = " << angle << endl; 
float EcalBarrel_depth=41.3;

float EcalBarrel_rmax= ( (EcalBarrel_rmin + EcalBarrel_depth + 1.0)/(TMath::Cos(3.14159/CaloSides))); 
float SolenoidBarrelInnerRadius1=EcalBarrel_rmax + 12.0;
//cout << " EcalBarrel_rmax = " << EcalBarrel_rmax << endl;

float SolenoidBarrelInnerCryostatThickness1=3.0;
float SolenoidBarrelInnerAirgapThickness1=1.5;
float SolenoidBarrelConductorInnerRadius1=SolenoidBarrelInnerRadius1 + 2.*SolenoidBarrelInnerCryostatThickness1 + SolenoidBarrelInnerAirgapThickness1;

float EcalPlug_rmin=Radius_SITLayer0+6.0;
  
float rmin=EcalPlug_rmin;
float rmax=SolenoidBarrelConductorInnerRadius1-10.;
float area= 3.14*(rmax*rmax-rmin*rmin);
int HcalPlug_fwd_layers=70+80+150;
float sumarea= (HcalPlug_fwd_layers*area)/1e4;//m^2
cout << " si area of hcalplug-fwd = " << sumarea << endl; 


//------------------HcalPlug_bwd-------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

float rmin1=EcalPlug_rmin;
float rmax1=SolenoidBarrelConductorInnerRadius1-10.;
int HcalPlug_bwd_layers=25+60+80;
float sumarea1=(3.14*(rmax1*rmax1-rmin1*rmin1)*HcalPlug_bwd_layers)/1e4;
cout << " si area of HcalPlug_bwd = " << sumarea1 << endl;


//-----------------EcalPlug_fwd------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------


float rmin2=EcalPlug_rmin; 
float rmax2=EcalBarrel_rmax;
//cout << " rmin2 " << rmin2 << endl; 
//cout << " rmax2 " << rmax2 << endl;
int  EcalPlug_fwd_layers=1+20+28;
float sumarea2=(3.14*(rmax2*rmax2-rmin2*rmin2)*EcalPlug_fwd_layers)/1e4;

cout << " si area of EcalPlug-fwd = " << sumarea2 << endl; 


//-------------------EcalPlug_bwd----------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

float rmin3=EcalPlug_rmin;
float rmax3=EcalBarrel_rmax;
int EcalPlug_bwd_layers=1+28+20;
float sumarea3=(3.14*(rmax3*rmax3-rmin3*rmin3)*EcalPlug_bwd_layers)/1e4;
cout << " si area of EcalPlug-bwd = " << sumarea3 << endl;

float sumall=sumarea+sumarea1+sumarea2+sumarea3;
cout << " -----------------total Si area for Lhe det. (m^2) = " << sumall << endl;//total Silicon area 
float number_of_reaout_channel=sumall*47059;
cout << " -----------------number of readout channels for LHe det. = " << number_of_reaout_channel << endl;//number of readout channel


/*******************************calculations of scintillator area**********************************/

//---------------------scintillator area of EcalBarrel---------------------------------------------
//-------------------------------------------------------------------------------------------------

float scint_area1;
float scint_area1t;
float scint_area2;
float scint_area2t;
float scint_area3;
float scint_area3t;

float Btd_disk_z4=-200.0; //cm
float Ftd_disk_z6 = 370.0;//cm
float EcalBarrel_zmax=Ftd_disk_z6-Btd_disk_z4+10;//!!

float Kapton_thickness= 0.03;//cm
float Polys_thickness=0.5;
float Pb_thickness1= 0.3;
float Air_thickness=0.15;

float ESpacing1 =(Polys_thickness+Air_thickness+Pb_thickness1);
float EfirstLength =EcalBarrel_rmin+Kapton_thickness+Polys_thickness;

float scint_area0=  (EcalBarrel_zmax*(2*3.14159*(EfirstLength)))/1e4;
cout << " scint_area0 for EcalBarrel = " << scint_area0 << endl;

int NbOfELayers1=9;

	for(int i=0; i!=NbOfELayers1; ++i)
	{
	scint_area1=(EcalBarrel_zmax*(2*3.14159*(EfirstLength+Pb_thickness1+Polys_thickness+i*ESpacing1)))/1e4;
	scint_area1t=scint_area1t+scint_area1;
         //cout << " scint_area1 for EcalBarrel = " << scint_area1t << endl;	
	}
cout << " scint_area1 for EcalBarrel = " << scint_area1t << endl;//m^2

int NbOfELayers2=14;
float Pb_thickness2= 0.4;

float FirstEcalvolume_length= NbOfELayers1*ESpacing1;
float ESpacing2 =(Polys_thickness+Air_thickness+Pb_thickness2);

	for(int i=0; i!=NbOfELayers2; ++i)
	{
	scint_area2=(EcalBarrel_zmax*(2*3.14159*(EfirstLength+FirstEcalvolume_length+Pb_thickness2+Polys_thickness+i*ESpacing2)))/1e4;
        scint_area2t=scint_area2t+scint_area2;
	//cout << " scint_area2 for EcalBarrel = " << scint_area2t << endl;	
	}
cout << " scint_area2 for EcalBarrel = " << scint_area2t << endl;//m^2

int NbOfELayers3=14;
float Pb_thickness3= 0.6;

float SecondEcalvolume_length= NbOfELayers2*ESpacing2;
float ESpacing3 =(Polys_thickness+Air_thickness+Pb_thickness3);

	for(int i=0; i!=NbOfELayers3; ++i)
	{
	scint_area3=(EcalBarrel_zmax*(2*3.14159*(EfirstLength+FirstEcalvolume_length+SecondEcalvolume_length+Pb_thickness3+Polys_thickness+i*ESpacing3)))/1e4;
        scint_area3t=scint_area3t+scint_area3;
	//cout << " scint_area2 for EcalBarrel = " << scint_area2 << endl;	
	}
cout << " scint_area3 for EcalBarrel = " << scint_area3t << endl;//m^2


float scint_area_for_EcalBarrel=scint_area0+scint_area1t+scint_area2t+scint_area3t;
cout << "scint_area_for_EcalBarrel(m^2)= " << scint_area_for_EcalBarrel << endl; //total scintillator area for EcalBarrel


//-----------------------------scintillator area of HcalBarrel-------------------------------------
//-------------------------------------------------------------------------------------------------

float HcalBarrel_length=EcalBarrel_zmax;

float scint_area_Hcal1;
float scint_area_Hcal1t;
float scint_area_Hcal2;
float scint_area_Hcal2t;
float scint_area_Hcal3;
float scint_area_Hcal3t;

float SolenoidBarrelQuenchbackThickness1=5.0;//cm
float SolenoidBarrelAlConductorThickness1=6.0;
float DipoleBarrelAlConductorThickness=1.5;
float SolenoidBarrelOuterAirgapThickness1 =3.0;
float SolenoidBarrelOuterCryostatThickness1=4.0;

float SolenoidBarrelOuterCryostatInnerRadius1=SolenoidBarrelConductorInnerRadius1  + SolenoidBarrelAlConductorThickness1 + SolenoidBarrelQuenchbackThickness1 + DipoleBarrelAlConductorThickness + SolenoidBarrelInnerCryostatThickness1 + 4.;

float SolenoidBarrelOuterRadius1=SolenoidBarrelOuterCryostatInnerRadius1 + SolenoidBarrelOuterAirgapThickness1 + SolenoidBarrelOuterCryostatThickness1;
float HcalBarrel_rmin=SolenoidBarrelOuterRadius1+5.0;//cm

int HcalBarrel_layers1=8;
int HcalBarrel_layers2=16;
int HcalBarrel_layers3=21;

float Steel235_thickness1=2.2;//cm
float Steel235_thickness2=2.5;//cm
float Steel235_thickness3=3.0;//cm

float Polys_thickness_H=0.50;
float Air_thickness_H=0.15;
float HSpacing1=Steel235_thickness1+Polys_thickness_H+Air_thickness_H;


for(int i=0; i!=HcalBarrel_layers1; ++i)
{scint_area_Hcal1=(HcalBarrel_length*(2*3.14159*(HcalBarrel_rmin+Steel235_thickness1+Polys_thickness_H+i*HSpacing1)))/1e4;
scint_area_Hcal1t=scint_area_Hcal1t+scint_area_Hcal1;
//cout << " scint_area_Hcal = " << scint_area_Hcalt << endl;	
}
cout << " scint_area_Hcal1 = " << scint_area_Hcal1t << endl;


float FirstHcalvolume_length= HcalBarrel_layers1*HSpacing1;
float HSpacing2 =Steel235_thickness2+Polys_thickness_H+Air_thickness_H;

for(int i=0; i!=HcalBarrel_layers2; ++i)
{scint_area_Hcal2=(HcalBarrel_length*(2*3.14159*(HcalBarrel_rmin+FirstHcalvolume_length+Steel235_thickness2+Polys_thickness_H+i*HSpacing2)))/1e4;
scint_area_Hcal2t=scint_area_Hcal2t+scint_area_Hcal2;
//cout << " scint_area_Hcal = " << scint_area_Hcalt << endl;	
}
cout << " scint_area_Hcal2 = " << scint_area_Hcal2t << endl;

float SecondHcalvolume_length= HcalBarrel_layers2*HSpacing2;
float HSpacing3 =Steel235_thickness3+Polys_thickness_H+Air_thickness_H;

for(int i=0; i!=HcalBarrel_layers3; ++i)
{scint_area_Hcal3=(HcalBarrel_length*(2*3.14159*(HcalBarrel_rmin+FirstHcalvolume_length+SecondHcalvolume_length+Steel235_thickness3+Polys_thickness_H+i*HSpacing3)))/1e4;
scint_area_Hcal3t=scint_area_Hcal3t+scint_area_Hcal3;
//cout << " scint_area_Hcal = " << scint_area_Hcalt << endl;	
}
cout << " scint_area_Hcal3 = " << scint_area_Hcal3t << endl;

float scint_area_for_HcalBarrel=scint_area_Hcal1t+scint_area_Hcal2t+scint_area_Hcal3t;
cout << "scint_area_for_HcalBarrel(m^2)= " << scint_area_for_HcalBarrel << endl; //total scintillator area for HcalBarrel


//-----------------------------scintillator area of HcalEndcap_fwd-------------------------------------
//-------------------------------------------------------------------------------------------------



float HcalBarrel_thickness=21.2+50.4+76.65;
float HcalEndcap_rmax = (HcalBarrel_rmin+HcalBarrel_thickness)/(TMath::Cos(3.14159/CaloSides));


float rmin4=HcalBarrel_rmin;
float rmax4=HcalEndcap_rmax;
int HcalEndcap_fwd_layers=10+20+28;
float sum_area1=(3.14*(rmax4*rmax4-rmin4*rmin4)*HcalEndcap_fwd_layers)/1e4;
cout << " scint. area of HcalEndcap_fwd = " << sum_area1 << endl;

//-----------------------------scintillator area of HcalEndcap_bwd-------------------------------------
//-------------------------------------------------------------------------------------------------



float rmin5=HcalBarrel_rmin;
float rmax5=HcalEndcap_rmax;
int HcalEndcap_bwd_layers=10+15+25;
float sum_area2=(3.14*(rmax5*rmax5-rmin5*rmin5)*HcalEndcap_bwd_layers)/1e4;
cout << " scint. area of HcalEndcap_bwd = " << sum_area2 << endl;









sumscint=scint_area_for_EcalBarrel+scint_area_for_HcalBarrel+sum_area1+sum_area2 ;
cout << " -----------------total sint. area for Lhe det. (m^2) = " << sumscint << endl;//total scintillator area for LHe det.
return 0;

}

