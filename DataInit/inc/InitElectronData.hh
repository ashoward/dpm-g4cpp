#ifndef InitElectronData_HH
#define InitElectronData_HH



class ElectronData;

void   InitElectronData(ElectronData& elData,
                        double electronCutInEnergy, double gammaCutInEnergy,
                        double parShigh, double parSlow, double parEcross);


void   InitElasticData(ElectronData& elData);
double FindScreeingParameter(double val);

void   InitElossData(ElectronData& elData, double electronCutInEnergy, double gammaCutInEnergy);

void   InitScatteringData(ElectronData& elData, double parShigh, double parSlow, double parEcross);

void   InitGSData(ElectronData& elData, double electronCutInEnergy);





#endif // InitElectronData_HH
