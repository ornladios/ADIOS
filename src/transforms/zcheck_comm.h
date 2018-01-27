static int zcheck_type(const int adios_type)
{
    int zc_type;
    switch (adios_type)
    {
    case adios_double:
        zc_type = ZC_DOUBLE;
        break;
    case adios_real:
        zc_type = ZC_FLOAT;
        break;
    default:
        adios_error(err_transform_failure, "No supported data type\n");
        return -1;
        break;
    }
    return zc_type;
}

static void zcheck_dim(
    const struct adios_file_struct *fd,
    const struct adios_var_struct *var,
    size_t r[5])
{
    struct adios_dimension_struct *d;
    d = var->pre_transform_dimensions;
    int ndims = (uint)count_dimensions(d);
    int i;
    for (i = 0; i < ndims; i++)
    {
        uint dsize = (uint)adios_get_dim_value(&d->dimension);
        if (fd->group->adios_host_language_fortran == adios_flag_yes)
            r[i] = dsize;
        else
            r[ndims - i - 1] = dsize;
        d = d->next;
    }
}

static void zcheck_write(
    ZC_DataProperty *dataProperty,
    ZC_CompareData *compareResult,
    struct adios_file_struct *fd,
    struct adios_var_struct *var)
{
#ifdef HAVE_ZCHECKER
    if (adios_verbose_level>7) ZC_printDataProperty(dataProperty);
    if (adios_verbose_level>7) ZC_printCompressionResult(compareResult);

    /*
    printf("psnr: %g\n", compareResult->psnr);
    printf("compressTime: %g\n", compareResult->compressTime);
    printf("compressRate: %g\n", compareResult->compressRate);
    printf("compressSize: %lu\n", compareResult->compressSize);
    printf("compressRatio: %g\n", compareResult->compressRatio);
    printf("decompressTime: %g\n", compareResult->decompressTime);
    printf("decompressRate: %g\n", compareResult->decompressRate);
    printf("minRelErr: %g\n", compareResult->minRelErr);
    printf("avgRelErr: %g\n", compareResult->avgRelErr);
    printf("maxRelErr: %g\n", compareResult->maxRelErr);
    printf("rmse: %g\n", compareResult->rmse);
    printf("nrmse: %g\n", compareResult->nrmse);
    printf("snr: %g\n", compareResult->snr);
    printf("pearsonCorr: %g\n", compareResult->pearsonCorr);
    */

    double my_entropy = dataProperty->entropy;
    double my_psnr = compareResult->psnr;
    double my_ratio = compareResult->compressRatio;
    double my_compressTime = compareResult->compressTime;
    double my_decompressTime = compareResult->decompressTime;
    log_debug("entropy, psnr, ratio: %g %g %g\n", my_entropy, my_psnr, my_ratio);
    int comm_size = 1;
    int comm_rank = 0;

    int64_t m_adios_file = (int64_t)fd;
    int64_t m_adios_group = (int64_t)fd->group;
    int64_t varid;

    char zname[255];
    sprintf(zname, "%s/%s", var->name, "entropy");
    varid = adios_common_define_var(m_adios_group, zname, "", adios_double, "", "", "");
    common_adios_write_byid(fd, (struct adios_var_struct *)varid, &my_entropy);

    sprintf(zname, "%s/%s", var->name, "psnr");
    varid = adios_common_define_var(m_adios_group, zname, "", adios_double, "", "", "");
    common_adios_write_byid(fd, (struct adios_var_struct *)varid, &my_psnr);

    sprintf(zname, "%s/%s", var->name, "ratio");
    varid = adios_common_define_var(m_adios_group, zname, "", adios_double, "", "", "");
    common_adios_write_byid(fd, (struct adios_var_struct *)varid, &my_ratio);

    sprintf(zname, "%s/%s", var->name, "compress_time");
    varid = adios_common_define_var(m_adios_group, zname, "", adios_double, "", "", "");
    common_adios_write_byid(fd, (struct adios_var_struct *)varid, &my_compressTime);

    sprintf(zname, "%s/%s", var->name, "decompress_time");
    varid = adios_common_define_var(m_adios_group, zname, "", adios_double, "", "", "");
    common_adios_write_byid(fd, (struct adios_var_struct *)varid, &my_decompressTime);
#else
    log_debug("%s: %s\n", "Z-checker", "Not available");
#endif
}
