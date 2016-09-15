#include <iostream>
#include <unittest++/UnitTest++.h>
#include <MeshTallyReader.hxx>
#include <McnpTally.hxx>

////////////////////////////
/// Helper Methods
////////////////////////////

const float MYPI = 3.1415;
const int MYDIM = 27;


McnpTally* MakeMcnpTally(float factor)
{
    McnpTally* aTally = new McnpTally;


    for(int i=0; i<MYDIM; i++){
        aTally->AddXBinTick(i);
        aTally->AddYBinTick(i*2);
        aTally->AddZBinTick(i*3);
    }

    aTally->AddEnergyTick(0.3);
    aTally->AddEnergyTick(1.2);
    aTally->InitValueArrays();

    aTally->SetTallyNumber(14);

    for(int i=0; i<MYDIM; i++)
        for(int j=0; j<MYDIM; j++)
            for(int k=0; k<MYDIM; k++)
                aTally->SetValue(0, i, j, k, factor*MYPI*float(i)*float(j)*float(k));

    return aTally;
}

McnpTally* MakeMcnpTally()
{
    return MakeMcnpTally(1.0f);
}

// /////////////////////////
///  Testing McnpTally
// /////////////////////////
TEST(Normalization)
{
    cout << "Testing normalization ...\n";
    // Create a McnpTally
    McnpTally* aTally = MakeMcnpTally();
    float factor(5.63);
    aTally->Multiply(factor);

    for(int i=0; i<MYDIM; i++)
        for(int j=0; j<MYDIM; j++)
            for(int k=0; k<MYDIM; k++){
                CHECK_CLOSE(MYPI*factor*float(i)*float(j)*float(k), aTally->GetValue(0,float(i),float(j),float(k)), 0.1f);
            }
}

TEST(SetAndGetTallyNumber)
{
    cout << "Testing Set and Get Tally Number ...\n";
    McnpTally* aTally = MakeMcnpTally();
    aTally->SetTallyNumber(99);
    CHECK_EQUAL(99, aTally->GetTallyNumber());
}

TEST(SetAndGetValue)
{
    cout << "Testing Set and Get Value ...\n";
    McnpTally* aTally = MakeMcnpTally();

    for(int i=0; i<MYDIM; i++)
        for(int j=0; j<MYDIM; j++)
            for(int k=0; k<MYDIM; k++)
            {
                aTally->SetValue(0, i, j, k, MYPI*float(i)*float(j)*float(k));
            }

    for(int i=0; i<MYDIM; i++)
        for(int j=0; j<MYDIM; j++)
            for(int k=0; k<MYDIM; k++)
            {
               // cout << i << " " << j << " " << k << " : " <<  MYPI*float(i)*float(j)*float(k) << " -- " << aTally->GetValue(0,i,j,k) << endl;
                CHECK_EQUAL(MYPI*float(i)*float(j)*float(k), aTally->GetValue(0,i,j,k));
            }
}


// //////////////////////////
///  Testing MeshTallyReader
// //////////////////////////
TEST(Operations)
{
    McnpTally* aTally = MakeMcnpTally();
    McnpTally* oTally = MakeMcnpTally(2.f);

    aTally->SetTallyNumber(1);
    oTally->SetTallyNumber(2);

    MeshTallyReader mtReader;
    mtReader.AddTally(*aTally);
    mtReader.AddTally(*oTally);

    string s_aMT("1.0");
    string s_oMT("2.0");
    string s_opAdd("+"), s_opMult("*"), s_opDiv("/");
    string array_Ops[3] = {s_opAdd, s_opMult, s_opDiv};

    int numNewTal(100);

    for(int o=0; o<3; o++)
    {
        cout << "Testing operation : " << array_Ops[o] << " ... " << endl;
        mtReader.PerformOperation(s_aMT, array_Ops[o], s_oMT);
        McnpTally newTally = mtReader.GetTally(numNewTal);

        CHECK_EQUAL(numNewTal, newTally.GetTallyNumber());
        numNewTal++;

        for(int i=0; i<newTally.XDimension(); i++)
            for(int j=0; j<newTally.YDimension(); j++)
                for(int k=0; k<newTally.ZDimension(); k++)
                {
                    if(o==0)
                         CHECK_CLOSE(aTally->GetValue(0,i,j,k)+oTally->GetValue(0,i,j,k), newTally.GetValue(0,i,j,k), 0.1);
                    else if(o==1)
                        CHECK_CLOSE(aTally->GetValue(0,i,j,k)*oTally->GetValue(0,i,j,k), newTally.GetValue(0,i,j,k), 0.1);
                   else if(o==2)
                    {
                        if(oTally->GetValue(0,i,j,k) == 0.f)
                        {
                            CHECK_CLOSE(0.f, newTally.GetValue(0,i,j,k), 0.1);
                        }
                        else
                        {
                            CHECK_CLOSE(aTally->GetValue(0,i,j,k)/oTally->GetValue(0,i,j,k), newTally.GetValue(0,i,j,k), 0.1);
                        }
                    }
                }
    }
}

// //////////////
/// Main
// ///////////////
int main(int argc, char* argv[])
{
    return UnitTest::RunAllTests();
}

