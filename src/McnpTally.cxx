#include "McnpTally.hxx"
#include <iostream>

McnpTally::McnpTally() 
{
    myArraysNotInitialized = true;
    myXSize = myYSize = myZSize = 1;
}


McnpTally::~McnpTally()
{
  /*delete[] myEnergyValueArray;
  myEnergyValueArray = NULL;
  delete[] myEnergyErrorArray;
  myEnergyValueArray = NULL;*/
}


void McnpTally::SetValue(int energyBin, int xBin, int yBin, int zBin, float value)
{
    if(myArraysNotInitialized)
    {
        std::cout << "ERROR : Arrays have not been initialized\n";
        return;
    }

    myEnergyValueArray[Position(energyBin, xBin, yBin, zBin)] = value;
}

float McnpTally::GetValue(int energyBin, int xBin, int yBin, int zBin) const
{
    return myEnergyValueArray[ Position(energyBin, xBin, yBin, zBin) ];
}


void McnpTally::SetError(int energyBin, int xBin, int yBin, int zBin, float error)
{
    if(myArraysNotInitialized)
        return;

    myEnergyErrorArray[Position(energyBin, xBin, yBin, zBin)] = error;
}

float McnpTally::GetError(int energyBin, int xBin, int yBin, int zBin) const
{
    return myEnergyErrorArray[Position(energyBin, xBin, yBin, zBin)];
}


void McnpTally::Multiply(float factor)
{
    for(int e=0; e<myESize; e++)
        for(int i=0; i<XDimension(); i++)
            for(int j=0; j<YDimension(); j++)
                for(int k=0; k<ZDimension(); k++ ){
                    myEnergyValueArray[Position(e , i , j , k) ] *= factor;
                }
}


void McnpTally::InitValueArrays()
{
    int energyBins = EnergyDimension();

    if(energyBins > 2) // additional array for SUM of all Energy Bins
    {
        energyBins++;
        float totalEnergy(-1.0f); // -1.0 is an avatar for the total energy range
        myEnergyBins.push_back(totalEnergy);
    }

    energyBins--;

    myESize = energyBins;
    myXSize = XDimension();
    myYSize = YDimension();
    myZSize = ZDimension();

    if( myESize*myXSize*myYSize*myZSize < 1 )
    {
        cout << "ERROR : Initialization of Arrays failed - array size equals ZERO!!!\n";
        return;
    }

    myEnergyValueArray = new float[myESize*myXSize*myYSize*myZSize];
    myEnergyErrorArray = new float[myESize*myXSize*myYSize*myZSize];

    // initialize
    for(int e=0; e< energyBins; e++)
    {
        for(int i=0; i< myXSize; i++)
        {
            for(int j=0; j< myYSize; j++)
            {
                for(int k=0; k< myZSize; k++)
                {
                    // cout << e << " " << i << " " << j << " " <<  k  << " ==> " << Position(e, i, j, k) << endl;
                    myEnergyValueArray[Position(e, i, j, k)] = -1.0f;
                    myEnergyErrorArray[Position(e, i, j, k)] = -1.0f;
                }
            }
        }
    }

    myArraysNotInitialized = false;
}


void McnpTally::SetMinVal(float energy, float val)
{
    list<float>::iterator itE = myEnergyBins.begin();
    itE++;
    unsigned int cnt(1);

    for(;itE != myEnergyBins.end(); itE++)
    {
        if((*itE) == energy)
            break;
        cnt++;
    }

    if(myLowerValueBoundary.size() < cnt)
        myLowerValueBoundary.push_back(val);
    else
    {
        list<float>::iterator itV = myLowerValueBoundary.begin();
        unsigned int cnt2(1);
        for(; itV != myLowerValueBoundary.end(); itV++)
        {
            if(cnt2 == cnt)
                break;
            cnt2++;
        }

        (*itV) = val;
    }
}

float McnpTally::GetMinVal(float energy) const
{
	if(myLowerValueBoundary.size() == 0)
		return 1.0e36f;

    list<float>::const_iterator itE = myEnergyBins.begin();
    itE++;
    list<float>::const_iterator itV = myLowerValueBoundary.begin();

    for(; itE != myEnergyBins.end(); itE++, itV++)
        if( *itE == energy )
            break;

    return *itV;
}

void McnpTally::SetMaxVal(float energy, float val)
{
    list<float>::iterator itE = myEnergyBins.begin();
    itE++;
    unsigned int cnt(1);

    for(;itE != myEnergyBins.end(); itE++)
    {
        if((*itE) == energy)
            break;
        cnt++;
    }

    if(myUpperValueBoundary.size() < cnt)
        myUpperValueBoundary.push_back(val);
    else
    {
        list<float>::iterator itV = myUpperValueBoundary.begin();
        unsigned int cnt2(1);
        for(; itV != myUpperValueBoundary.end(); itV++)
        {
            if(cnt2 == cnt)
                break;
            cnt2++;
        }

        (*itV) = val;
    }
}

float McnpTally::GetMaxVal(float energy) const
{
	if(myUpperValueBoundary.size() == 0)
		return 0.0f;

    list<float>::const_iterator itE = myEnergyBins.begin();
    itE++;
    list<float>::const_iterator itV = myUpperValueBoundary.begin();

    for(; itE != myEnergyBins.end(); itE++, itV++)
        if( *itE == energy )
            break;

    return *itV;
}


float McnpTally::GetEnergyTick(int number) const
{
    list<float>::const_iterator itE = myEnergyBins.begin();
//    itE++;
    int cnt(0);

    for(; itE != myEnergyBins.end(); itE++, cnt++)
        if(cnt == number)
            break;

    return *itE;

}


int McnpTally::Position(int e, int i, int j, int k) const
{
    return (k + j*myZSize + i*myZSize*myYSize + e*myZSize*myYSize*myXSize);
}
