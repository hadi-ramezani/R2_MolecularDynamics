#ifndef BONDED_H
#define	BONDED_H

class Initial;
class Parameters;
class Vector;
struct BondElem;
struct AngleElem;
struct DihedralElem;
struct ImproperElem;


class Bonded {
public:
    // constructor creates the lists of bonds, angles, etc.
    Bonded(const Initial *, Parameters *);
    ~Bonded();

    // compute the energy, given the set of coordinates
    void compute(const Vector *coords, Vector *f, double& Ebond, double& Eangle,
                   double &Edihedral, double &Eimproper) const;

private:
    int nbonds;
    int nangles;
    int ndihedrals;
    int nimpropers;

    BondElem *bonds;
    AngleElem *angles; 
    DihedralElem *dihedrals;
    ImproperElem *impropers;
   
    void build_bondlist(const Initial *, const Parameters *);
    void build_anglelist(const Initial *, const Parameters *);
    void build_dihedrallist(const Initial *, const Parameters *);
    void build_improperlist(const Initial *, const Parameters *);

    double compute_bonds(const Vector *, Vector *) const;
    double compute_angles(const Vector *, Vector *) const;
    double compute_dihedrals(const Vector *, Vector *) const;
    double compute_impropers(const Vector *, Vector *) const;
};

#endif	/* BONDED_H */

