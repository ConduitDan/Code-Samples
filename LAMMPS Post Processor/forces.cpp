#include "forces.h"
double Geometry::area(std::vector<LammpsAtom> vertices, std::vector<std::array<int, 3>> facets) {
    std::array<double,3> v1;
    std::array<double,3> v2;
    std::array<double,3> crossProduct;
    double area = 0;
    #pragma omp parallel for private(v1,v2,crossProduct) reduction(+:area)
    for (auto face = facets.begin(); face<facets.end(); ++face){
        v1 = minus(vertices[(*face)[1]].get_position(),vertices[(*face)[0]].get_position());
        v2 = minus(vertices[(*face)[2]].get_position(),vertices[(*face)[1]].get_position());
        crossProduct = cross(v1,v2);
        area +=sqrt(dot(crossProduct,crossProduct))/2;
    }

    return area;
}
std::array<double,3> Geometry::cross(std::array<double,3> v1, std::array<double,3> v2){
// | i    j    k    |
// | v1_0 v1_1 v1_2 |
// | v2_0 v2_1 v2_2 |
    std::array<double,3> out;
    out[0] = v1[1]*v2[2]-v1[2]*v2[1];
    out[1] = -v1[0]*v2[2]+v2[0]*v1[2];
    out[2] = v1[0]*v2[1]-v1[1]*v2[0];
    return out;
}
std::array<double,3> Geometry::minus(std::array<double,3> v1, std::array<double,3> v2){
    std::array<double,3> out;
    out[0] = v1[0]-v2[0];
    out[1] = v1[1]-v2[1];
    out[2] = v1[2]-v2[2];
    return out;
}
double Geometry::dot(std::array<double,3> v1, std::array<double,3> v2){
    double out = v1[0]*v2[0];
    out+=v1[1]*v2[1];
    out+=v1[2]*v2[2];
    return out;
}
int Geometry::bond_diff(std::vector<std::array<int, 3>> bonds1, std::vector<std::array<int, 3>> bonds2) {
    // just go though and compare 1 to one, if this doesn't work we have to copy them to sets and try to set compare
    int sum = 0;
    typedef std::unordered_set<int> inner_type;
    typedef std::unordered_set<inner_type, hash_on_set_of_ints> set_of_unique_sets;


    set_of_unique_sets bonds1Set;
    set_of_unique_sets bonds2Set;
    for (int i = 0; i<bonds1.size(); i++){
        bonds1Set.insert({bonds1[i][0],bonds1[i][1]});
        bonds2Set.insert({bonds2[i][0],bonds2[i][1]});
    }
    
    for (int i = 0; i<bonds1.size(); i++){
        if (bonds1[i]!=bonds2[i]) {
            sum++;
        }
    }
    return sum;
}
double Geometry::length(LammpsAtom v1, LammpsAtom v2) {
    double dx = v1.get_position()[0]-v2.get_position()[0];
    double dy = v1.get_position()[1]-v2.get_position()[1];
    double dz = v1.get_position()[2]-v2.get_position()[2];
    return sqrt(dx*dx+dy*dy+dz*dz); 
}
double Geometry::radius(std::vector<LammpsAtom> vertices, int type) {
    double r = 0;
    int count = 0;
    #pragma omp parallel for reduction(+:r,count)
    for (auto vertex = vertices.begin(); vertex<vertices.end(); ++vertex){
        if (vertex->get_type()!=type) continue;
        r += sqrt(dot(vertex->get_position(),vertex->get_position()));
        count++;
    }
    return r/count;
}

double Potentials::sw_harmonic(sw_harmonic_params params, std::array<int, 2> bond, std::vector<LammpsAtom>* atoms) {
    std::array<double,3> r0 = (*atoms)[bond[0]-1].get_position();
    std::array<double,3> r1 = (*atoms)[bond[1]-1].get_position();

    std::array<double,3> r10 = Geometry::minus(r1,r0);
    double rsq = Geometry::dot(r10,r10);
    double r = sqrt(rsq);

    double E = 0;
    if (r<params.rmin||r>params.rmax) {
        return std::numeric_limits<double>::infinity();
    }
    if (r<params.rmin+params.delta) E+= params.sigma*exp(1.0/(r-(params.rmin+params.delta)))/(r-params.rmin);
    if (r>params.rmax-params.delta) E+= params.sigma*exp(1.0/((params.rmax-params.delta)-r))/(params.rmax-r);
    E+=(r-params.leq)*(r-params.leq)*0.5*params.k;
    return E;

}
double Potentials::dihedral(double stiffness, std::array<int, 4> dihedral, std::vector<LammpsAtom>* atoms) {
  int a = dihedral[0]-1;
  int b = dihedral[1]-1;
  int c = dihedral[2]-1;
  int d = dihedral[3]-1;

  double xa = (*atoms)[a].get_position()[0];
  double ya = (*atoms)[a].get_position()[1];
  double za = (*atoms)[a].get_position()[2];

  double xb = (*atoms)[b].get_position()[0];
  double yb = (*atoms)[b].get_position()[1];
  double zb = (*atoms)[b].get_position()[2];

  double xc = (*atoms)[c].get_position()[0];
  double yc = (*atoms)[c].get_position()[1];
  double zc = (*atoms)[c].get_position()[2];

  double xd = (*atoms)[d].get_position()[0];
  double yd = (*atoms)[d].get_position()[1];
  double zd = (*atoms)[d].get_position()[2];

  // Calculate angle theta which is the angle bw the two planes made by
  // {a,c,b} and {b,c,d}

  double dx_cb = xc - xb;
  double dy_cb = yc - yb;
  double dz_cb = zc - zb;

  double dx_ab = xa - xb;
  double dy_ab = ya - yb;
  double dz_ab = za - zb;

  double dx_dc = xd - xc;
  double dy_dc = yd - yc;
  double dz_dc = zd - zc;

  double N1_x = (dy_cb * dz_dc) - (dz_cb * dy_dc);
  double N1_y = (dz_cb * dx_dc) - (dx_cb * dz_dc);
  double N1_z = (dx_cb * dy_dc) - (dy_cb * dx_dc);

  double N2_x = (-dy_cb * dz_ab) - (-dz_cb * dy_ab);
  double N2_y = (-dz_cb * dx_ab) - (-dx_cb * dz_ab);
  double N2_z = (-dx_cb * dy_ab) - (-dy_cb * dx_ab);

  double norm_N1 = std::sqrt(N1_x * N1_x + N1_y * N1_y + N1_z * N1_z);
  double norm_N2 = std::sqrt(N2_x * N2_x + N2_y * N2_y + N2_z * N2_z);

  N1_x = N1_x / norm_N1;
  N1_y = N1_y / norm_N1;
  N1_z = N1_z / norm_N1;

  N2_x = N2_x / norm_N2;
  N2_y = N2_y / norm_N2;
  N2_z = N2_z / norm_N2;

  double costheta = N1_x * N2_x + N1_y * N2_y + N1_z * N2_z;
  
  
  // Returns the energy associated with the dihedral
  return stiffness*(1-costheta);

}

double Potentials::power(power_params params, double l){

    if (l<params.lmin){
        return params.k*std::pow(params.lmin-l,params.r-1)*std::pow(params.r,2+params.r)/l;
    }
    if (l>params.lmax){
        return -params.k*std::pow(l-params.lmax,params.r-1)*std::pow(params.r,2+params.r)/l;
    }
    return 0;
}
double Potentials::power_harmonic(power_params params, double l){
    double f = 0;
    if (l<params.lmin){
        f += params.k*std::pow(params.lmin-l,params.r-1)*std::pow(params.r,2+params.r)/l;
    }
    if (l>params.lmax){
        f+= -params.k*std::pow(l-params.lmax,params.r-1)*std::pow(params.r,2+params.r)/l;
    }
    f += -params.k*(l-(params.lmax+params.lmin)/2)/l;
    return f;
}
