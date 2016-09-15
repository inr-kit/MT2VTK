#ifndef MCNPMESHGEOMETRYTYPE_HXX
#define MCNPMESHGEOMETRYTYPE_HXX

enum McnpMeshFormat{
    Mcnp_IJ,
    Mcnp_IK,
    Mcnp_JK,
    Mcnp_COL,
    Mcnp_COL_WITH_ENERGY,
    Mcnp_CF
};

enum McnpMeshGeometry{
    Mcnp_CARTESIAN,
    Mcnp_CYLINDICAL
};

#endif // MCNPMESHGEOMETRYTYPE_HXX
