#ifndef TURQMD_HD
#define TURQMD_HD
#include <string>
#include <cstring>
#include <stdexcept>
extern "C" {
#include <stdbool.h>
}
#include <cstdio>
#include <iostream>

namespace TURQMD{

const int nmax = 40000; // maximum number of particles
const int smax = 500; // maximum number of spectators
const int numcto = 400; // maximum number of CTOption
const int numctp = 400; // maximum number of CTParam
extern "C" {
    void c_urqmd_main_();
    void c_urqmd_init_();
    void c_urqmd_fin_();
    void c_urqmd_evt_(int*);

    void c_set_pro_(int*, int*);
    void c_set_pro2_(int*, int*);
    void c_set_tar_(int*, int*);
    void c_set_tar2_(int*, int*);
    void c_set_box_(double*, double*, int*, int*);
    void c_set_bpet_(int*, int*, int*, int*, double*);
    void c_set_eneelb_(int*, double*);
    void c_set_plb_(double*);
    void c_set_plbplg_(int*, double*, double*, int*);
    void c_set_ecm_(double*);
    void c_set_eneelg_(int*, double*, double*, int*);
    void c_set_imp_(double*);
    void c_set_imp2_(double*, double*);
    void c_set_eos_(int*);
    void c_set_nev_(int*);
    void c_set_rsd_(int*);
    void c_set_cdt_(double*);
    void c_set_tim_(double*, double*);
    void c_set_stb_(int*);
    void c_supf13_();
    void c_supf14_();
    void c_supf15_();
    void c_supf16_();
    void c_supf19_();
    void c_supf20_();
    void c_use_inputfile_();
    void c_old_event_();
    void c_use_file17_();

    void c_on_external_seed_();
    void c_off_external_seed_();
    void c_set_evt_seed_(int*);
    void c_get_evt_seed_(int*);
    void c_s_cto_(int*, double*);
    void c_g_cto_(int*, double*);
    void c_s_ctp_(int*, double*);
    void c_g_ctp_(int*, double*);
    void c_s_nevents_(double*);
    void c_g_nevents_(double*);
    void c_s_npart_(int*);
    void c_g_npart_(int*);
    void c_s_cnpartit_(int*);
    void c_g_cnpartit_(int*);

    void c_wdata_(int*, float*,float*,float*);

    extern struct { // particles' coordination common block
        double r0[nmax], rx[nmax], ry[nmax], rz[nmax],
               p0[nmax], px[nmax], py[nmax], pz[nmax],
               fmass[nmax], rww[nmax], dectime[nmax];
    } coor_;

    extern struct {
        int spin[nmax],ncoll[nmax],charge[nmax],ityp[nmax],
            lstcoll[nmax],iso3[nmax],origin[nmax],uid[nmax];
    } isys_;

    extern struct { // spectators' common block
        double r0s[smax], rxs[smax], rys[smax], rzs[smax],
               p0s[smax], pxs [smax],pys[smax], pzs[smax],
               sfmass[smax];
    } scoor_;

    extern struct {
        int spin[smax], scharge[smax], sityp[smax],
            siso3[smax], suid[smax];
    } sisys_;

    extern struct {
        // number of spectators in an event
        int nspec;
    } ssys_;

    extern struct { // options' common block
        int CTOption[numcto], CTParam[numctp];
    } options_;

    extern struct {
      bool fixedseed,bf13,bf14,bf15,bf16,bf17,
           bf19, bf20;
    } loptions_;

    extern struct { // urqmd system's common block
        int npart, nbar, nmes, ctag, nsteps,uid_cnt,
            ranseed,event,Ap,At,Zp,Zt,eos,dectag,
            NHardRes,NSoftRes,NDecRes,NElColl,NBlColl;
        bool success;
    } sys_;

    extern struct {
        double time,acttime,bdist,bimp,bmin,ebeam,ecm;
    } rsys_;

    extern struct { // input common block
        int nevents,spityp,prspflg,trspflg,
            spiso3,outsteps,bflag,srtflag,efuncflag,nsrt,
            firstev, npb;
    } inputs_;

    extern struct {
        double srtmin,srtmax,pbeam,betann,betatar,betapro,
               pbmin, pbmax;
    } input2_;

    extern struct {
        int cnpartit, cmul;
        int c_current_is_empty;
    } curqmd_rout_;

    extern struct {
        int c_quiet, c_skip_empty_event;
    } curqmd_opt_;
} // end of extern 'C' for fortran lib

///// Setting of Text Option
int Ap, Zp, At, Zt, spityp_1, spityp_2, spiso3_1, spiso3_2;
int opt;
double lbox,edens;
int solid,para,mbox;
double bdist, bmin;
double ebeam, pbeam;
double pmin, pmax;
int npb;
double Ecm;
int nsrt;
double srtmin, srtmax, dtimestep, caltim, outtim;
int efuncflag;
int EoS, nevents, ranseed, partid;
int bptpart_, bptityp_, bptiso3_;
double bptpmax_;
void setpro(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d %d", &Ap, &Zp);
    c_set_pro_(&Ap, &Zp);
}
void setpro2(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d %d", &spityp_1, &spiso3_1);
    c_set_pro2_(&spityp_1, &spiso3_1);
}
void settar(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d %d", &At, &Zt);
    c_set_tar_(&At, &Zt);
}
void settar2(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d %d", &spityp_2, &spiso3_2);
    c_set_tar2_(&spityp_2, &spiso3_2);
}
void setbox(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg %lg %d %d", &lbox, &edens, &solid, &para);
    c_set_box_(&lbox, &edens, &solid, &para);
}
void setbpe(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d %d %d %lg", &bptityp_,&bptiso3_,&bptpart_,&bptpmax_);
    opt = 1;
    c_set_bpet_(&opt, &bptityp_,&bptiso3_,&bptpart_,&bptpmax_);
}
void setbpt(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d %d %d", &bptityp_, &bptiso3_, &bptpart_);
    opt = 2;
    c_set_bpet_(&opt, &bptityp_, &bptiso3_, &bptpart_, &bptpmax_);
}
void setene(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg", &ebeam);
    opt = 1;
    c_set_eneelb_(&opt, &ebeam);
}
void setelb(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg", &ebeam);
    opt = 2;
    c_set_eneelb_(&opt, &ebeam);
}
void setplb(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg", &pbeam);
    c_set_plb_(&pbeam);
}
void setplb2(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg %lg %d", &pmin, &pmax, &npb);
    opt = 1;
    c_set_plbplg_(&opt, &pmin, &pmax, &npb);
}
void setplg(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg %lg %d", &pmin, &pmax, &npb);
    opt = 2;
    c_set_plbplg_(&opt, &pmin, &pmax, &npb);
}
void setecm(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg", &Ecm);
    c_set_ecm_(&Ecm);
}
void setene2(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg %lg %d", &srtmin, &srtmax, &nsrt);
    opt = 1;
    c_set_eneelg_(&opt, &srtmin, &srtmax, &nsrt);
}
void setelg(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg %lg %d", &srtmin, &srtmax, &nsrt);
    opt = 2;
    c_set_eneelg_(&opt, &srtmin, &srtmax, &nsrt);
}
void setimp(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg", &bdist);
    c_set_imp_(&bdist);
}
void setimp2(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg %lg", &bmin, &bdist);
    c_set_imp2_(&bmin, &bdist);
}
void seteos(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d", &EoS);
    c_set_eos_(&EoS);
}
void setnev(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d", &nevents);
    c_set_nev_(&nevents);
}
void setrsd(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d", &ranseed);
    c_set_rsd_(&ranseed);
}
void setcdt(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg", &dtimestep);
    c_set_cdt_(&dtimestep);
}
void settim(const std::string &sopt) {
    sscanf(sopt.c_str(), "%lg %lg", &caltim, &outtim);
    c_set_tim_(&caltim, &outtim);
}
void setstb(const std::string &sopt) {
    sscanf(sopt.c_str(), "%d", &partid);
    c_set_stb_(&partid);
}
///// End of Setting of Text Option

///// Class ValProxy
typedef void (*ValDeg0)(double*);
typedef void (*ValDeg1)(int*, double*);
typedef void (*ValDeg2)(int*, int*, double*);
typedef void (*iValDeg0)(int*);
typedef void (*iValDeg1)(int*, int*);
typedef void (*iValDeg2)(int*, int*, int*);

// TODO using template instead of overload
class ValProxy {
public:
    typedef void (ValProxy::*DegFunc)(double&) const;
    typedef void (ValProxy::*iDegFunc)(int&) const;
    // for double function
    ValProxy(ValDeg0 sHandle, ValDeg0 gHandle) {
        sdeg0 = sHandle;
        gdeg0 = gHandle;
        sdeg = &ValProxy::_sdeg_0d;
        gdeg = &ValProxy::_gdeg_0d;
        isdeg = &ValProxy::_sdeg_0d_int;
        igdeg = &ValProxy::_gdeg_0d_int;
    }
    ValProxy(ValDeg1 sHandle, ValDeg1 gHandle) {
        sdeg1 = sHandle;
        gdeg1 = gHandle;
        sdeg = &ValProxy::_sdeg_1d;
        gdeg = &ValProxy::_gdeg_1d;
        isdeg = &ValProxy::_sdeg_1d_int;
        igdeg = &ValProxy::_gdeg_1d_int;
    }
    ValProxy(ValDeg2 sHandle, ValDeg2 gHandle) {
        sdeg2 = sHandle;
        gdeg2 = gHandle;
        sdeg = &ValProxy::_sdeg_2d;
        gdeg = &ValProxy::_gdeg_2d;
        isdeg = &ValProxy::_sdeg_2d_int;
        igdeg = &ValProxy::_gdeg_2d_int;
    }
    // for int function
    ValProxy(iValDeg0 sHandle, iValDeg0 gHandle) {
        isdeg0 = sHandle;
        igdeg0 = gHandle;
        isdeg = &ValProxy::_isdeg_0d;
        igdeg = &ValProxy::_igdeg_0d;
        sdeg = &ValProxy::_isdeg_0d_double;
        gdeg = &ValProxy::_igdeg_0d_double;
    }
    ValProxy(iValDeg1 sHandle, iValDeg1 gHandle) {
        isdeg1 = sHandle;
        igdeg1 = gHandle;
        isdeg = &ValProxy::_isdeg_1d;
        igdeg = &ValProxy::_igdeg_1d;
        sdeg = &ValProxy::_isdeg_1d_double;
        gdeg = &ValProxy::_igdeg_1d_double;
    }
    ValProxy(iValDeg2 sHandle, iValDeg2 gHandle) {
        isdeg2 = sHandle;
        igdeg2 = gHandle;
        isdeg = &ValProxy::_isdeg_2d;
        igdeg = &ValProxy::_igdeg_2d;
        sdeg = &ValProxy::_isdeg_2d_double;
        gdeg = &ValProxy::_igdeg_2d_double;
    }
    const ValProxy& operator = (double val) const {
        (this->*sdeg)(val);
        return *this;
    }
    const ValProxy& operator = (int val) const {
        (this->*isdeg)(val);
        return *this;
    }
    const ValProxy& operator() (int idx) const {
        idx1 = idx;
        return *this;
    }
    const ValProxy& operator() (int idx, int idy) const {
        idx1 = idx; idx2 = idy;
        return *this;
    }
    operator double() const {
        double ret;
        (this->*gdeg)(ret);
        return ret;
    }
    int toint() const {
        int ret;
        (this->*igdeg)(ret);
        return ret;
    }
private:
    DegFunc sdeg, gdeg;
    iDegFunc isdeg, igdeg;
    ValDeg0 sdeg0, gdeg0;
    ValDeg1 sdeg1, gdeg1;
    ValDeg2 sdeg2, gdeg2;
    iValDeg0 isdeg0, igdeg0;
    iValDeg1 isdeg1, igdeg1;
    iValDeg2 isdeg2, igdeg2;
    mutable int idx1, idx2, idx3;
    // for double
    void _sdeg_0d(double &val) const {
        double v = val;
        sdeg0(&v);
    }
    void _sdeg_1d(double &val) const {
        int idx = idx1;
        double v = val;
        sdeg1(&idx, &v);
    }
    void _sdeg_2d(double &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        double v = val;
        sdeg2(&idx1, &idx2, &v);
    }
    void _gdeg_0d(double &val) const {
        double ret;
        gdeg0(&ret);
        val = ret;
    }
    void _gdeg_1d(double &val) const {
        int idx = idx1;
        double ret;
        gdeg1(&idx, &ret);
        val = ret;
    }
    void _gdeg_2d(double &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        double ret;
        gdeg2(&idx1, &idx2, &ret);
        val = ret;
    }
    ///
    void _sdeg_0d_int(int &val) const {
        double v = (int) val;
        sdeg0(&v);
    }
    void _sdeg_1d_int(int &val) const {
        int idx = idx1;
        double v = (int) val;
        sdeg1(&idx, &v);
    }
    void _sdeg_2d_int(int &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        double v = (int) val;
        sdeg2(&idx1, &idx2, &v);
    }
    void _gdeg_0d_int(int &val) const {
        double ret;
        gdeg0(&ret);
        val = (int) ret;
    }
    void _gdeg_1d_int(int &val) const {
        int idx = idx1;
        double ret;
        gdeg1(&idx, &ret);
        val = (int) ret;
    }
    void _gdeg_2d_int(int &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        double ret;
        gdeg2(&idx1, &idx2, &ret);
        val = (int) ret;
    }
    // for int
    void _isdeg_0d(int &val) const {
        int v = val;
        isdeg0(&v);
    }
    void _isdeg_1d(int &val) const {
        int idx = idx1;
        int v = val;
        isdeg1(&idx, &v);
    }
    void _isdeg_2d(int &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        int v = val;
        isdeg2(&idx1, &idx2, &v);
    }
    void _igdeg_0d(int &val) const {
        int ret;
        igdeg0(&ret);
        val = ret;
    }
    void _igdeg_1d(int &val) const {
        int idx = idx1;
        int ret;
        igdeg1(&idx, &ret);
        val = ret;
    }
    void _igdeg_2d(int &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        int ret;
        igdeg2(&idx1, &idx2, &ret);
        val = ret;
    }
    ///
    void _isdeg_0d_double(double &val) const {
        int v = (int) val;
        isdeg0(&v);
    }
    void _isdeg_1d_double(double &val) const {
        int idx = idx1;
        int v = (int) val;
        isdeg1(&idx, &v);
    }
    void _isdeg_2d_double(double &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        int v = (int) val;
        isdeg2(&idx1, &idx2, &v);
    }
    void _igdeg_0d_double(double &val) const {
        int ret;
        igdeg0(&ret);
        val = (double) ret;
    }
    void _igdeg_1d_double(double &val) const {
        int idx = idx1;
        int ret;
        igdeg1(&idx, &ret);
        val = (double) ret;
    }
    void _igdeg_2d_double(double &val) const {
        int idx1 = idx1;
        int idx2 = idx2;
        int ret;
        igdeg2(&idx1, &idx2, &ret);
        val = (double) ret;
    }
};
///// End of Class ValProxy

////// Class StrOptProxy
typedef void (*StrOptDeg)(const std::string &sopt);

class StrOptProxy {
public:
    StrOptProxy(StrOptDeg sdegf) : strdegf(sdegf) {};
    std::string operator << (const std::string &sopt);
    std::string operator >> (std::string &sout);
private:
    StrOptDeg strdegf;
    std::string stropt;
};

std::string StrOptProxy::operator << (const std::string &sopt) {
    strdegf(sopt);
    stropt = sopt;
    return sopt;
}

std::string StrOptProxy::operator >> (std::string &sout) {
    sout = stropt;
    return sout;
}
////// End of Class StrOptProxy

///// Class TUrQMD
class TUrQMD;

typedef void (*ctimeFuncT)(TUrQMD&, float);
TUrQMD* _TUrQMD_tcallback_instance;
typedef void (*ValDeg)(int*, double*);

class TUrQMD  {
public:
    TUrQMD() :
        record_time_evolution(false),
        record_time_function(NULL),
        CTO(c_s_cto_, c_g_cto_),
        CTP(c_s_ctp_, c_g_ctp_),
        nevents(c_s_nevents_, c_g_nevents_),
        Npart(c_s_cnpartit_, c_g_cnpartit_),
        pro(setpro), PRO(setpro2), tar(settar),
        TAR(settar2), nev(setnev), tim(settim),
        ene(setene), elb(setelb),
        plb(setplb), PLB(setplb2),
        PLG(setplg), ecm(setecm),
        ENE(setene2), ELG(setelg),
        imp(setimp), IMP(setimp2),
        eos(seteos), box(setbox),
        bpt(setbpt), bpe(setbpe), rsd(setrsd),
        stb(setstb), cdt(setcdt),
        Npart_(curqmd_rout_.cnpartit),
        mul(curqmd_rout_.cmul),
        bimp(rsys_.bimp),
        quiet_out(curqmd_opt_.c_quiet),
        skip_empty_event(curqmd_opt_.c_skip_empty_event),
        is_empty_event(curqmd_rout_.c_current_is_empty)
    {
        c_urqmd_main_();
        use_external_seed = true;
        /*
        get_para();
        c_off_timeana_();
        _TJam_tcallback_instance = NULL;
        */
        // init text output file
        for (int i = 0; i < bfsize; ++i) {
            f[i] = false;
        }
    }
    void init() {
        //set_para();
        /*
        if (record_time_evolution) {
            c_on_timeana_();
            _TJam_tcallback_instance=this;
        }
        */
        //get_para();
        if (use_external_seed) {
            c_on_external_seed_();
        } else {
            c_off_external_seed_();
        }
        parse_input();
        c_urqmd_init_();
        _iev = 0;
    }
    void finish() {
        //set_para();
        c_urqmd_fin_();
        //get_para();
    }
    void evt(const int &i) {
        if (use_external_seed) {
            throw std::runtime_error("Please specify the random seed by using TUrQMD::evt(int eventCount, int seed)");
        }
        int iv = i;
        int cur_sed;
        //set_para();
        c_urqmd_evt_(&iv);
        if (skip_empty_event) {
            if (is_empty_event) {
                puts("retry event");
                evt(i);
            }
        }
        //get_para();
        c_get_evt_seed_(&cur_sed);
        current_event_seed = cur_sed;
        _iev = i;
    }
    void evt(const int &i, const int &nseed) {
        if (!use_external_seed) {
            throw std::runtime_error("Cannot specify the random seed event by event when TUrQMD::use_external_seed==false");
        }
        int iv = i;
        int seed = nseed, cur_sed;
        //set_para();
        c_set_evt_seed_(&seed);
        c_urqmd_evt_(&iv);
        if (skip_empty_event) {
            if (is_empty_event) {
                puts("retry event");
                //TODO nseed maybe duplicated
                evt(i, nseed+1);
            }
        }
        //get_para();
        c_get_evt_seed_(&cur_sed);
        current_event_seed = cur_sed;
        _iev = i;
    }
    void write_test() {
        int i;
        using namespace std;
        puts("CTO:");
        for (i = 0; i < 10; ++i) {
            cout << options_.CTOption[i] << " ";
        }
        puts("");
        puts("rx:");
        for (i = 0; i < 10; ++i) {
            cout << coor_.rx[i] << " ";
        }
        puts("");
    }
    void write_pdata(int *pid, float *px, float *py, float *pz) {
        c_wdata_(pid, px, py, pz);
    }
    /*
    void write_pxdata(int *Npart, int *ks, int *pid, float *rx, float *ry, float *rz, float *rt, float *rft,
            float *px, float *py, float *pz, float *e, float *mass, float *charge) {
        c_wxdata_(Npart, ks, pid, rx, ry, rz, rt, rft, px, py, pz, e, mass, charge);
    }
    void write_time(float* catime, int* Npart, int* cmul, int *ks_ary, int* pid_ary,
            float* rx_ary, float* ry_ary, float* rz_ary, float *rt_ary, float *rft_ary,
            float* px_ary, float* py_ary, float* pz_ary,
            float* e_ary, float* m_ary, float* chge_ary)
        c_tphase_(catime, Npart, cmul, ks_ary, pid_ary,
            rx_ary, ry_ary, rz_ary, rt_ary, rft_ary,
            px_ary, py_ary, pz_ary,
            e_ary, m_ary, chge_ary);
    }
    */
    static const int bfsize = 50;
    int iev() const { return _iev; }
    int maxv() const { return 30000; }
    int mevent, nstep, mxcoll;
    bool f[bfsize];
    double bmin, bmax, dt, nv, nbary, nmeson, mentry;
    std::string cwin, frame, proj, targ;
    StrOptProxy pro, PRO, tar, TAR, nev, tim;
    StrOptProxy ene, elb, plb, PLB, PLG, ecm;
    StrOptProxy ENE, ELG, imp, IMP, eos, box;
    StrOptProxy bpt, bpe, rsd, stb, cdt;
    bool record_time_evolution;
    bool use_external_seed;
    ctimeFuncT record_time_function;
    ValProxy CTO, CTP, nevents, Npart;
    int &Npart_;
    int &mul;
    double &bimp;
    int &quiet_out, &is_empty_event, &skip_empty_event;
private:
    char _cwin[20], _frame[20], _proj[20], _targ[20];
    int _iev, current_event_seed;
    void parse_input() {
        if (f[9]) {
            c_use_inputfile_();
        }
        if (f[10]) {
            c_old_event_();
        }
        if (!f[13]) {
            c_supf13_();
        }
        if (!f[14]) {
            c_supf14_();
        }
        if (!f[15]) {
            c_supf15_();
        }
        if (!f[16]) {
            c_supf16_();
        }
        if (f[17]) {
            c_use_file17_();
        }
        if (!f[19]) {
            c_supf19_();
        }
        if (!f[20]) {
            c_supf20_();
        }
    }
    /*
    void set_para() {
        cstr_ck();
        c_s_para_(&mevent, &bmin, &bmax, _cwin, _frame, _proj, _targ, &dt, &nstep);
    }
    void get_para() {
        c_para_(&mevent, &bmin, &bmax, _cwin, _frame, _proj, _targ, &dt, &nstep, &nv, &nbary, &nmeson, &mxcoll, &mentry);
        cstr_rt();
    }
    */
};
///// End of Class TUrQMD

/*
extern "C" {
void c_jam_time_record_(double *atime) {
    _TJam_tcallback_instance->record_time_function(*_TJam_tcallback_instance, float(*atime));
}
}
*/

} // end of namespace TURQMD

#endif
