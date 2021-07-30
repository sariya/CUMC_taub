#ifndef FUNC_create_regenie_H
#define FUNC_create_regenie_H
//Date 07 30 2021

using namespace std;
#include <string>

struct snp_store
{
    std::string snp_name; //store chr12:1234:REF_allele:ALT_allele
    std::string func_annotation; // store exoninc, intronic and such infor
} ;

struct gene
{
    std::string gene_name;  //abca7
    std::string chromosome; // 9
    std::string start_position; // 123
};

void remove_trailingspaces(char *newline);

#endif
// typedef struct snp
// {
//     std::string snp_name;
//     std::string func_annotation;
// } snp_store;

// typedef struct gene_SNP
// {
//     std::string gene_name;
// } gene_snp_store;
