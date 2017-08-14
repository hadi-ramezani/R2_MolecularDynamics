#include <math.h>
#include <iostream>

#include "LJTable.h"
#include "Parameters.h"

LJTable::LJTable(const Parameters *params)
{
  table_dim = params->get_num_vdw_params();
  half_table_sz = table_dim * table_dim;

  table = new LJTableEntry[half_table_sz * 2];

  for (int i=0; i < table_dim; i++)
    for (int j=i; j < table_dim; j++)
    {
      LJTableEntry *curij = &(table[i*table_dim+j]);
      LJTableEntry *curji = &(table[j*table_dim+i]);
      compute_vdw_params(params,i,j,curij,curij+half_table_sz);

      // Copy to transpose entry
      *curji = *curij;
      *(curji + half_table_sz) = *(curij + half_table_sz);
    }
}

LJTable::~LJTable()
{
  delete [] table;
}

void LJTable::compute_vdw_params(const Parameters *params, int i, int j,
				 LJTableEntry *cur, 
				 LJTableEntry *cur_scaled) {
  double A, B, A14, B14;
  double sigma_max;
  //  We need the A and B parameters for the Van der Waals.  These can
  //  be explicitly be specified for this pair or calculated from the
  //  sigma and epsilon values for the two atom types
  if (params->get_vdw_pair_params(i,j, &A, &B, &A14, &B14))
  {
    cur->A = A;
    cur->B = B;
    cur_scaled->A = A14;
    cur_scaled->B = B14;

    double sigma_ij = pow(A/B,1./6.);
    double sigma_ij14 = pow(A14/B14,1./6.);

    sigma_max = ( sigma_ij > sigma_ij14 ? sigma_ij : sigma_ij14 );
  }
  else
  {
    //  We didn't find explicit parameters for this pair. So instead,
    //  get the parameters for each atom type separately and use them
    //  to calculate the values we need
    double sigma_i, sigma_i14, epsilon_i, epsilon_i14;
    double sigma_j, sigma_j14, epsilon_j, epsilon_j14;

    params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
				       &epsilon_i14,i);
    params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14, 
				       &epsilon_j14,j);
  	
    double sigma_ij = 0.5 * (sigma_i+sigma_j);
    double epsilon_ij = sqrt(epsilon_i*epsilon_j);
    double sigma_ij14 = 0.5 * (sigma_i14+sigma_j14);
    double epsilon_ij14 = sqrt(epsilon_i14*epsilon_j14);

    sigma_max = ( sigma_ij > sigma_ij14 ? sigma_ij : sigma_ij14 );

    //  Calculate sigma^6
    sigma_ij *= sigma_ij*sigma_ij;
    sigma_ij *= sigma_ij;
    sigma_ij14 *= sigma_ij14*sigma_ij14;
    sigma_ij14 *= sigma_ij14;
    
    //  Calculate LJ constants A & B
    cur->B = 4.0 * sigma_ij * epsilon_ij;
    cur->A = cur->B * sigma_ij;
    cur_scaled->B = 4.0 * sigma_ij14 * epsilon_ij14;
    cur_scaled->A = cur_scaled->B * sigma_ij14;
  }
  //  Calculate exclcut2
  cur_scaled->exclcut2 = cur->exclcut2 = 0.64 * sigma_max * sigma_max;

}

