I=phantom(100);
writecfl('/autofs/space/daisy_002/users/Gilad/A',I);
B=I*0+2;
C=I*0+3;
D=I*0+4;
E=I*0+5;
writecfl('/autofs/space/daisy_002/users/Gilad/B',B);
writecfl('/autofs/space/daisy_002/users/Gilad/C',C);
writecfl('/autofs/space/daisy_002/users/Gilad/D',D);
writecfl('/autofs/space/daisy_002/users/Gilad/E',E);
%%
bart('linopScript')

X=readcfl('/autofs/space/daisy_002/users/Gilad/X');
grmss(X(abs(I)>0)./I(abs(I)>0))

Y=readcfl('/autofs/space/daisy_002/users/Gilad/Y');
grmss(Y(abs(I)>0)./I(abs(I)>0))
%%
debug_printf(DP_INFO,"Test\n");
    long ADims[DIMS];
    complex float* A = load_cfl("/autofs/space/daisy_002/users/Gilad/A", DIMS, ADims);
    printf("Reading %ld %ld %ld\n",ADims[0],ADims[1],ADims[2]);

    complex float* B = load_cfl("/autofs/space/daisy_002/users/Gilad/B", DIMS, ADims);
    complex float* C = load_cfl("/autofs/space/daisy_002/users/Gilad/C", DIMS, ADims);
    complex float* D = load_cfl("/autofs/space/daisy_002/users/Gilad/D", DIMS, ADims);
    complex float* E = load_cfl("/autofs/space/daisy_002/users/Gilad/E", DIMS, ADims);
    debug_printf(DP_INFO,"Read ABCD\n");
    complex float* X = create_cfl("/autofs/space/daisy_002/users/Gilad/X", DIMS, ADims);
    complex float* Y = create_cfl("/autofs/space/daisy_002/users/Gilad/Y", DIMS, ADims);

    long AFlags=md_nontriv_dims(DIMS,ADims);
    
    const struct linop_s* fop1 = linop_fmac_create(DIMS, ADims, ~AFlags, ~AFlags, ~AFlags, B);
    const struct linop_s* fop2 = linop_fmac_create(DIMS, ADims, ~AFlags, ~AFlags, ~AFlags, C);
    struct fmac_data * pdata=getpdata();

    const struct linop_s* fop3 = linop_fmac_create(DIMS, ADims, ~AFlags, ~AFlags, ~AFlags, D);

    const struct linop_s* fop=linop_chain(fop1,fop2);
    fop=linop_chain(fop,fop3);

    linop_forward(fop, DIMS, ADims, X,	DIMS, ADims, A);

    debug_printf(DP_INFO,"OK X\n");
    setfmacdataToNewTensor(pdata, E);

    linop_forward(fop, DIMS, ADims, Y,	DIMS, ADims, A);
    debug_printf(DP_INFO,"OK Y\n");



    linop_free(fop1);
    linop_free(fop2);
    linop_free(fop3);
    linop_free(fop);
    debug_printf(DP_INFO,"Freed ops\n");

    unmap_cfl(DIMS, ADims, X);
    unmap_cfl(DIMS, ADims, Y);
    unmap_cfl(DIMS, ADims, A);
    unmap_cfl(DIMS, ADims, B);
    unmap_cfl(DIMS, ADims, C);
    unmap_cfl(DIMS, ADims, D);
    unmap_cfl(DIMS, ADims, E);
    debug_printf(DP_INFO,"Freed ABCD,X\n");
    debug_printf(DP_INFO,"end test\n");
    debug_printf(DP_INFO,"end test\n");
    debug_printf(DP_INFO,"end test\n");
    debug_printf(DP_INFO,"end test\n");