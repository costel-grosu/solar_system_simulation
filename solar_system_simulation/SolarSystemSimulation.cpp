

#include "stdafx.h"

#include "Platform.h"
#include "vector3.h"

#include "BasicIntegrators.h"


static const double G = 0.01720209895; // KG
static const double C = 173.144632720536344565; // C  Speed of light, au / d
static const double AU = 1.4959787066e11; // Astronomical unit in meters

typedef enum { Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Sun } CelestialBody;

class NBody {
public:
    int n;
    V3 *x, *v, *a, *jerk;
    V3 *xp, *vp, *ap, *jp;
    double *m, *m_GR;

    bool useRelativity;

    NBody(int n, bool useRelativity) {
        this->n = n;
        x = new V3[n];
        v = new V3[n];
        a = new V3[n];
        jerk = new V3[n];
        xp = new V3[n];
        vp = new V3[n];
        ap = new V3[n];
        jp = new V3[n];
        m = new double[n];
        m_GR = new double[n];
        this->useRelativity = useRelativity;
    }

    ~NBody() {
        delete[] x;
        delete[] v;
        delete[] a;
        delete[] jerk;
        delete[] xp;
        delete[] vp;
        delete[] ap;
        delete[] jp;
        delete[] m;
        delete[] m_GR;
    }

    static NBody* initSolarSystem(bool useRelativity) {
        CelestialBody i;
        
        //init with JPL data
        //2440400.50 JD0 Julian date of initial state vector
        const FLOAT T0 = 2440400.50;
        
        NBody* ss = new NBody(10, useRelativity);

        i = Sun;
        ss->x[i] = V3(4.5144118714356666407e-003, 7.2282841152065867346e-004, 2.4659100492567986271e-004);
        ss->v[i] = V3(-2.8369446340813151639e-007, 5.1811944086463255444e-006, 2.2306588340621263489e-006);
        ss->m[i] = 2.95912208285591102582E-4;

        i = Mercury;
        ss->x[i] = V3(3.6030663368975339466e-001, -9.4812876741771684223e-002, -8.7466840117233140436e-002);
        ss->v[i] = V3(3.7085069798210382186e-003, 2.4854958767430945324e-002, 1.2929109014677844626e-002);
        ss->m[i] = 4.91254745145081175785E-11;

        i = Venus;
        ss->x[i] = V3(6.0786466491731583464e-001, -3.5518362463675619232e-001, -1.9824142909855515515e-001);
        ss->v[i] = V3(1.1156645711264016669e-002, 1.5494075513638794325e-002, 6.2773904546696609267e-003);
        ss->m[i] = 7.24345620963276523095E-10;

        i = Earth;
        ss->x[i] = V3(1.0820747754938311664e-001, -9.2711110739430602933e-001, -4.0209347855944090112e-001);
        ss->v[i] = V3(1.6833251020051496668e-002, 1.5602036176919105255e-003, 6.7646174015847273137e-004);
        ss->m[i] = 8.88769273403302327042E-10;

        i = Mars;
        ss->x[i] = V3(-1.2796408611369531836e-001, -1.3262618005333617013e+000, -6.0530808652523961512e-001);
        ss->v[i] = V3(1.4481919298277924969e-002, 8.0528538390447499843e-005, -3.5188931029397090065e-004);
        ss->m[i] = 9.54952894222405763492E-11;

        i = Jupiter;
        ss->x[i] = V3(-5.3896824544609061333e+000, -7.7026549518616593034e-001, -1.9866431165907522014e-001);
        ss->v[i] = V3(1.0053452569924098185e-003, -6.5298425191689416643e-003, -2.8258787532429609536e-003);
        ss->m[i] = 2.82534210344592625472E-7;

        i = Saturn;
        ss->x[i] = V3(7.9527768530257360864e+000, 4.5078822184006553686e+000, 1.5201955253183338898e+000);
        ss->v[i] = V3(-3.1594662504012930114e-003, 4.3714634278372622354e-003, 1.9441395169137103763e-003);
        ss->m[i] = 8.45946850483065929285E-8;

        i = Uranus;
        ss->x[i] = V3(-1.8278236586353147533e+001, -9.5764572881482056433e-001, -1.6132190397271035415e-001);
        ss->v[i] = V3(1.7108310564806817248e-004, -3.7646704682815900043e-003, -1.6519678610257000136e-003);
        ss->m[i] = 1.28881623813803488851E-8;

        i = Neptune;
        ss->x[i] = V3(-1.6367191358770888335e+001, -2.3760896725373076342e+001, -9.3213866179497290101e+000);
        ss->v[i] = V3(2.6225242764289213785e-003, -1.5277473123858904045e-003, -6.9183197562182804864e-004);
        ss->m[i] = 1.53211248128427618918E-8;

        i = Pluto;
        ss->x[i] = V3(-3.0447680255169362534e+001, -5.3177934960261367037e-001, 9.0596584886274922101e+000);
        ss->v[i] = V3(2.8177758090360373050e-004, -3.1469590804946202045e-003, -1.0794238049289112837e-003);
        ss->m[i] = 2.27624775186369921644E-12;

        return ss;
    }

    // given x, v calculate a, jerk
    void updateForces(V3* x, V3* v, V3* a, V3* jerk) {
        int i, j;
        V3 dx, dv, da, dj;
        double r2, ir2, ir3;
        for (i = 0; i < n; i++) {
            a[i] = V3(0, 0, 0);
            jerk[i] = V3(0, 0, 0);
        }

        for (i = 0; i < n; i++)
            for (j = i+1; j < n; j++) {
                dx = x[i] - x[j];
                dv = v[i] - v[j];
                ir2 = 1 / (dx * dx);
                ir3 = ir2 * sqrt(ir2);
                
                dj = (dv + dx * (-3 * (dx * dv) * ir2));

                a[i] -= dx * (m[j] * ir3); 
                a[j] += dx * (m[i] * ir3);
                jerk[i] -= dj * (m[j] * ir3);
                jerk[j] += dj * (m[i] * ir3);
            }

        if (useRelativity) {
            relativity(x, v, a);
        }
    }

    // update acceleration with relativity correction
    void relativity(V3* x, V3* v, V3* a) {
        int i, j, k;
        V3 t0, t1, t2, t3, t4;
        double invc2 = 1 / (C * C);

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) 
                if (i != j) {
                    V3 r_ij = x[i] - x[j];
                    double invr2_ij = 1 / (r_ij * r_ij);
                    double invr_ij = sqrt(invr2_ij);
                    V3 n_ij = r_ij * invr_ij;

                    // post-newtonian formula from "Gravity  Newtonian, Post-newtonian, Relativistic" page 441
                    V3 t0 = -m[j] * invr2_ij * n_ij; // newtonian acceleration
                    t1 = t0 * (v[i]*v[i] - 4.f*(v[i]*v[j]) + 2.f*(v[j]*v[j]) - 3.f/2.f*(n_ij*v[j]) - 5.f*m[i]*invr_ij - 4.f*m[j]*invr_ij);
                    t2 = -t0 * (n_ij*(4.f*v[i] - 3.f*v[j])) * (v[i] - v[j]);

                    // t3,t4 seem to be 100 times smaller than t1, t2
                    t3 = V3(0, 0, 0);
                    t4 = V3(0, 0, 0);
                    for (k = 0; k < n; k++)
                        if (k != i && k != j) {
                            V3 r_jk = x[j] - x[k];
                            double invr2_jk = 1 / (r_jk * r_jk);
                            double invr_jk = sqrt(invr2_jk);
                            V3 n_jk = r_jk * invr_jk;

                            t3 += m[j] * m[k] * invr2_ij * (4.f*invr_ij + invr_jk - invr2_jk / (2.f*invr_ij)*(n_ij*n_jk)) * n_ij;
                            t4 -= 7.f / 2.f * m[j] * m[k] * invr2_jk * invr_ij * n_jk;
                        }
                    //a[i] += invc2 * (t1 + t2);
                    a[i] += invc2 * (t1 + t2 + t3 + t4);
                }
        }
    }

    // given x, v calculate acc, jerk
    void updateForcesFast(V3* x, V3* v, V3* a, V3* jerk) {
        int i, j;
        double dx, dy, dz, dvx, dvy, dvz, dax, day, daz, djx, djy, djz;
        double r2, ir2, ir3;

        for (i = 0; i < n; i++) {
            a[i] = V3(0, 0, 0);
            jerk[i] = V3(0, 0, 0);
        }

        for (i = 0; i < n; i++)
            for (j = i + 1; j < n; j++) {
                dx = x[i].x - x[j].x;
                dy = x[i].y - x[j].y;
                dz = x[i].z - x[j].z;
                dvx = v[i].x - v[j].x;
                dvy = v[i].y - v[j].y;
                dvz = v[i].z - v[j].z;

                ir2 = 1 / (dx * dx + dy * dy + dz * dz);
                ir3 = ir2 * sqrt(ir2);
                
                double dxdv = dx * dvx + dy *dvy + dz * dvz;

                djx = (dvx + dx * (-3 * dxdv * ir2));
                djy = (dvy + dy * (-3 * dxdv * ir2));
                djz = (dvz + dz * (-3 * dxdv * ir2));

                double mir3 = m[i] * ir3; 
                double mjr3 = m[j] * ir3;

                a[i].x -= dx * mjr3;
                a[i].y -= dy * mjr3;
                a[i].z -= dz * mjr3;

                a[j].x += dx * mir3;
                a[j].y += dy * mir3;
                a[j].z += dz * mir3;

                jerk[i].x -= djx * mjr3;
                jerk[i].y -= djy * mjr3;
                jerk[i].z -= djz * mjr3;

                jerk[j].x += djx * mir3;
                jerk[j].y += djy * mir3;
                jerk[j].z += djz * mir3;
            }

        if (useRelativity) {
            relativity(x, v, a);
        }
    }

    //4th order hermite scheme
    void integrate(double dt) {
        int i, j;
        double dt2 = dt * dt;

        updateForces(x, v, a, jerk);

        //predictor
        for (i = 0; i < n; i++) {
            xp[i] = x[i] + v[i] * dt + a[i] * (dt2 / 2) + jerk[i] * (dt2 * dt / 6);
            vp[i] = v[i] + a[i] * dt + jerk[i] * (dt2 / 2);
        }
        updateForces(xp, vp, ap, jp);

        // corrector
        for (i = 0; i < n; i++) {
            V3 oldv = v[i];
            v[i] = v[i] + (a[i] + ap[i]) * (dt / 2) + (jerk[i] - jp[i]) * (dt2 / 12);
            x[i] = x[i] + (v[i] + oldv) * (dt / 2) + (a[i] - ap[i]) * (dt2 / 12);
        }
    }


    double energy() {
        // calculate potential + kinetic energy
        int i, j;
        double K = 0; // kinetic
        double U = 0; // potential
        for (i = 0; i < n; i++) {
            double const Ki = (v[i] * v[i]) * (m[i] / 2);
            K += Ki;
            for (j = i+1; j < n; j++) {
                V3 r = x[i] - x[j];
                double const Uij = (m[i] * m[j]) / sqrt(r * r);
                U -= Uij;
            }
        }
        // remove the extra G (each mass is Gm)
        return(K + U)/G;
    }


    // constrain the barycenter to stay at the origin, move the sun.
    void barycenterFixsun() {
        int i, j, k;
        double ic2 = 1 / (C * C);
        for (k = 0; k < 2; k++) {
            for (i = 0; i < n; i++) {
                double s = 0;
                for (j = i+1; j < n; j++) {
                    V3 dx = x[i] - x[j];
                    s -= m[j] / sqrt(dx * dx);
                }
                m_GR[i] = m[i] * (1 + ((v[i] * v[i]) * 0.5 + s) * ic2);
            }

            V3 cx = V3(0, 0, 0);
            V3 cv = V3(0, 0, 0);
            for (i = 0; i < n; i++) {
                if (i != Sun) {
                    cx += x[i] * m_GR[i];
                    cv += v[i] * m_GR[i];
                }
            }
            printf("barycenter move %g\n", x[Sun] + cx / m_GR[Sun]);
            x[Sun] = -cx / m_GR[Sun];
            v[Sun] = -cv / m_GR[Sun];
        }
    }

    static void testSolarSystem() {
        NBody * ss = NBody::initSolarSystem(true);
        double dt;
        double t;
#if 0
        { //for perf measurement
            __int64 t0 = now();
            for (int i = 0; i < 1000; i++) {
                s.updateForces(s.x, s.v, s.a, s.jerk);
            }
            __int64 t1 = now();
            printf("ticks for updateForce: %d\n", (int)((t1 - t0) / 1000));
        }

        { //for perf measurement
            __int64 t0 = now();
            for (int i = 0; i < 1000; i++) {
                s.updateForcesFast(s.x, s.v, s.a, s.jerk);
            }
            __int64 t1 = now();
            printf("ticks for updateForceFast: %d\n", (int)((t1 - t0) / 1000));
        }
#endif

        double e0 = ss->energy();
        printf("E(t0)=%g\n", e0);
        dt = 0.1;
        int nsteps = 0;
        //ss->barycenterFixsun();
        __int64 t0 = now();
        for (t = 2440400.50; t <= 2440800.50; t += dt) {
            nsteps++;
            ss->integrate(dt);
        }
        __int64 t1 = now();
        //s->barycenterFixsun();

        printf("at %g after %d steps, relative delta Energy=%g\n", 
            t, nsteps, (ss->energy()-e0)/e0);

        printf("ticks for loop: %d\n", (int)((t1 - t0) / nsteps));

        { // check positional error
            V3 xtest, vtest, dx, dv;
            // test data at 2440800.50

            // some error for Mercury, need general relativity
            xtest = V3(-3.4350297408563828055e-001, -2.5118910315168408963e-001, -9.9170891227838064646e-002);
            vtest = V3(1.1640059959070938016e-002, -1.7964875644231269124e-002, -1.0817311005210173814e-002);
            dx = ss->x[Mercury] - xtest;
            printf("mercury dx=(%g %g %g)\n", dx.x, dx.y, dx.z);
            dv = ss->v[Mercury] - vtest;
            printf("mercury dv=(%g %g %g)\n", dv.x, dv.y, dv.z);

            // smaller error for Jupiter
            xtest = V3(-4.2373615722355993645e+000, -3.1464733091277224306e+000, -1.2462241659031373286e+000);
            vtest = V3(4.6228653666202124358e-003, -5.0555918692615562352e-003, -2.2818455722493645242e-003);
            dx = ss->x[Jupiter] - xtest;
            printf("jupiter dx=(%g %g %g)\n", dx.x, dx.y, dx.z);
            dv = ss->v[Jupiter] - vtest;
            printf("Jupiter dv=(%g %g %g)\n", dv.x, dv.y, dv.z);
        }
    }

};




int _tmain(int argc, _TCHAR* argv[])
{

    NBody::testSolarSystem();
//    Integrators::test();

    printf("Press a key to close.");
    getchar();

	return 0;
}

