#ifndef LJTABLE_H
#define LJTABLE_H

class Parameters;

struct LJTableEntry
{
  double exclcut2;
  double A;
  double B;
};

class LJTable
{
public:

  LJTable(const Parameters *);
  ~LJTable();

  const LJTableEntry *table_val_scaled14(int i, int j) const {
    return table + i * table_dim + j + half_table_sz;
  }
    
  const LJTableEntry *table_val(int i, int j) const {
    return table + i * table_dim + j;
  } 

private:

  void compute_vdw_params(const Parameters *, int i, int j, 
			  LJTableEntry *cur, LJTableEntry *cur_scaled);
  LJTableEntry *table;
  int half_table_sz;
  int table_dim;
};

#endif

