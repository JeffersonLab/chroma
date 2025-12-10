#define QLA_c_eq_r_div_c(r,a,b) { QLA_Real ab = a/QLA_norm2_c(b); QLA_c_eq_r_plus_ir(r, ab*QLA_real(b), -ab*QLA_imag(b)); }

// these need to have nv, lat and sub defined

#define create_V(r) { for(int _i=0; _i<nv; _i++) (r)[_i] = QDPN(create_V_L)(nc, lat); }
#define destroy_V(r) { for(int _i=0; _i<nv; _i++) QDPN(destroy_V)((r)[_i]); }

#define V_eq_zero(r) for(int _i=0; _i<nv; _i++) QDPN(V_eq_zero)((r)[_i],sub)

  //#define V_eq_V(r,a) QDPN(V_veq_V)(r,a,sub,nv)
#define V_eq_V(r,a) for(int _i=0; _i<nv; _i++) QDPN(V_eq_V)((r)[_i],(a)[_i],sub)

#define V_peq_V(r,a) for(int _i=0; _i<nv; _i++) QDPN(V_peq_V)((r)[_i],(a)[_i],sub)
#define V_meq_V(r,a) for(int _i=0; _i<nv; _i++) QDPN(V_meq_V)((r)[_i],(a)[_i],sub)

  //#define r_eq_norm2_V(r,a) { QLA_Real _r[nv]; QDPN(r_veq_norm2_V)(_r,a,sub,nv); *(r) = 0; for(int _i=0; _i<nv; _i++) {printf0("%g\n",_r[_i]); *(r) += _r[_i];} }
#define r_eq_norm2_V(r,a) { QLA_Real _r; *(r) = 0; for(int _i=0; _i<nv; _i++) { QDPN(r_eq_norm2_V)(&_r,(a)[_i],sub); *(r) += _r; } }

#define r_eq_re_V_dot_V(r,a,b) { QLA_Real _r; *(r) = 0; for(int _i=0; _i<nv; _i++) { QDPN(r_eq_re_V_dot_V)(&_r,(a)[_i],(b)[_i],sub); *(r) += _r; } }

    //#define c_eq_V_dot_V(r,a,b) { QLA_Complex _r[nv]; QDPN(c_veq_V_dot_V)(_r,a,b,sub,nv); QLA_c_eq_r(*r,0); for(int _i=0; _i<nv; _i++) QLA_c_peq_c(*r, _r[_i]); }
#define c_eq_V_dot_V(r,a,b) { QLA_Complex _r; QLA_c_eq_r(*(r),0); for(int _i=0; _i<nv; _i++) { QDPN(c_eq_V_dot_V)(&_r,(a)[_i],(b)[_i],sub); QLA_c_peq_c(*(r), _r); } }

#define V_eq_r_times_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_eq_r_times_V)((r)[_i],a,(b)[_i],sub); }
#define V_peq_r_times_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_peq_r_times_V)((r)[_i],a,(b)[_i],sub); }
#define V_meq_r_times_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_meq_r_times_V)((r)[_i],a,(b)[_i],sub); }

    //#define V_eq_c_times_V(r,a,b) { QLA_Complex _a[nv]; for(int _i=0; _i<nv; _i++) QLA_c_eq_c(_a[_i],*a); QDPN(V_veq_c_times_V)(r,_a,b,sub,nv); }
#define V_eq_c_times_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_eq_c_times_V)((r)[_i],a,(b)[_i],sub); }

    //#define V_peq_c_times_V(r,a,b) { QLA_Complex _a[nv]; for(int _i=0; _i<nv; _i++) QLA_c_eq_c(_a[_i],*a); QDPN(V_vpeq_c_times_V)(r,_a,b,sub,nv); }
#define V_peq_c_times_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_peq_c_times_V)((r)[_i],a,(b)[_i],sub); }

    //#define V_meq_c_times_V(r,a,b) { QLA_Complex _a[nv]; for(int _i=0; _i<nv; _i++) QLA_c_eq_c(_a[_i],*a); QDPN(V_vmeq_c_times_V)(r,_a,b,sub,nv); }
#define V_meq_c_times_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_meq_c_times_V)((r)[_i],a,(b)[_i],sub); }

#define V_eq_V_plus_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_eq_V_plus_V)((r)[_i],(a)[_i],(b)[_i],sub); }

#define V_eq_V_minus_V(r,a,b) { for(int _i=0; _i<nv; _i++) QDPN(V_eq_V_minus_V)((r)[_i],(a)[_i],(b)[_i],sub); }
