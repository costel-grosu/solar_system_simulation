// solar_system_simulation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include  <assert.h>
#include <cmath>
#include "vector3.h"

#define now __rdtsc

#define NOINLINE __declspec(noinline)
#define INLINE __declspec(inline)

#define FLOAT double
#define V3 dvec3

#define G 1.0

struct Body {
    V3 x, v, a;
    V3 jerk;

    FLOAT Gm;
    FLOAT t;

    // Predicted values
    V3 xp, vp, ap;
    V3 jerkp;
    FLOAT tp;

    // acceleration of m1 towards mg
    // mg already multiplied with G
    static V3 acceleration(V3 r, V3 v, FLOAT Gm1, FLOAT Gm2) {
        // F = G * m1 * m2 / r^2 * (R/r)
        // a = F / m1  -- so Newtonian acceleration does not depend on m1
        FLOAT r2 = r * r;
        FLOAT rabs = 1 / sqrt(r2);
        V3 A = (Gm2 / r2) * (r * rabs);
        return A;
    };


    // Hermite 4th order predictor
    void Predictor4(FLOAT dt) {
        const double dt2 = (1. / 2.) * dt;
        const double dt3 = (1. / 3.) * dt;
        xp = x + dt * (v + dt2 * (a + dt3 * jerk));
        vp = v + dt * (a + dt2 * jerk);
    };


    void Corrector4(FLOAT dt) {
        const double h = dt / 2;
        const double h2 = - h / 3;
        V3 vc = v + h*((ap + a) + h2 * (jerkp - jerk));
        x = x + h*((vc + v) + h2 * (ap - a));
        a = ap;
        jerk = jerkp;
        t = tp;
    }


};


class SolarSystem {
public:
    static const int NBody = 10;
    Body planets[NBody];

    void init(void) {
        //init with JPL data
        //2440400.50 JD0 Julian date of initial state vector
        const FLOAT T0 = 2440400.50;

         /* Sun */
        planets[0].x = V3(4.5144118714356666407e-003, 7.2282841152065867346e-004, 2.4659100492567986271e-004);
        planets[0].v = V3(-2.8369446340813151639e-007, 5.1811944086463255444e-006, 2.2306588340621263489e-006);
        planets[0].Gm = 2.95912208285591102582E-4;
        planets[0].t = 2440400.50;

        /* Mercury */
        planets[1].x = V3(3.6030663368975339466e-001, -9.4812876741771684223e-002, -8.7466840117233140436e-002);
        planets[1].v = V3(3.7085069798210382186e-003, 2.4854958767430945324e-002, 1.2929109014677844626e-002);
        planets[1].Gm = 4.91254745145081175785E-11;
        planets[1].t = 2440400.50;

        /* Venus */
        planets[2].x = V3(6.0786466491731583464e-001, -3.5518362463675619232e-001, -1.9824142909855515515e-001);
        planets[2].v = V3(1.1156645711264016669e-002, 1.5494075513638794325e-002, 6.2773904546696609267e-003);
        planets[2].Gm = 7.24345620963276523095E-10;
        planets[2].t = 2440400.50;

        /* EMB */
        planets[3].x = V3(1.0820747754938311664e-001, -9.2711110739430602933e-001, -4.0209347855944090112e-001);
        planets[3].v = V3(1.6833251020051496668e-002, 1.5602036176919105255e-003, 6.7646174015847273137e-004);
        planets[3].Gm = 8.88769273403302327042E-10;
        planets[3].t = 2440400.50;

        /* Mars */
        planets[4].x = V3(-1.2796408611369531836e-001, -1.3262618005333617013e+000, -6.0530808652523961512e-001);
        planets[4].v = V3(1.4481919298277924969e-002, 8.0528538390447499843e-005, -3.5188931029397090065e-004);
        planets[4].Gm = 9.54952894222405763492E-11;
        planets[4].t = 2440400.50;

        /* Jupiter */
        planets[5].x = V3(-5.3896824544609061333e+000, -7.7026549518616593034e-001, -1.9866431165907522014e-001);
        planets[5].v = V3(1.0053452569924098185e-003, -6.5298425191689416643e-003, -2.8258787532429609536e-003);
        planets[5].Gm = 2.82534210344592625472E-7;
        planets[5].t = 2440400.50;

        /* Saturn */
        planets[6].x = V3(7.9527768530257360864e+000, 4.5078822184006553686e+000, 1.5201955253183338898e+000);
        planets[6].v = V3(-3.1594662504012930114e-003, 4.3714634278372622354e-003, 1.9441395169137103763e-003);
        planets[6].Gm = 8.45946850483065929285E-8;
        planets[6].t = 2440400.50;

        /* Uranus */
        planets[7].x = V3(-1.8278236586353147533e+001, -9.5764572881482056433e-001, -1.6132190397271035415e-001);
        planets[7].v = V3(1.7108310564806817248e-004, -3.7646704682815900043e-003, -1.6519678610257000136e-003);
        planets[7].Gm = 1.28881623813803488851E-8;
        planets[7].t = 2440400.50;

        /* Neptune */
        planets[8].x = V3(-1.6367191358770888335e+001, -2.3760896725373076342e+001, -9.3213866179497290101e+000);
        planets[8].v = V3(2.6225242764289213785e-003, -1.5277473123858904045e-003, -6.9183197562182804864e-004);
        planets[8].Gm = 1.53211248128427618918E-8;
        planets[8].t = 2440400.50;

        /* Pluto */
        planets[9].x = V3(-3.0447680255169362534e+001, -5.3177934960261367037e-001, 9.0596584886274922101e+000);
        planets[9].v = V3(2.8177758090360373050e-004, -3.1469590804946202045e-003, -1.0794238049289112837e-003);
        planets[9].Gm = 2.27624775186369921644E-12;
        planets[9].t = 2440400.50;

        /* MOON */
        //planets[0].x = V3(-8.35703164195952601e-4, -1.98543915768166071e-3, -1.08326877048661754e-3);
        //planets[0].v = V3(5.98752118417335956e-4, -1.74153713527242722e-4, -8.84771962437116281e-5);
        // 1.09318924524284471369E-11
    };


    FLOAT energy() {
        // calculate potential + kinetic energy
        int i, j;
        FLOAT K = 0; // kinetic
        FLOAT U = 0; // potential
        for (i = 0; i < NBody; i++) {
            const Body *pi = &planets[i];
            FLOAT const Ki = (pi->v.norm2() * pi->Gm) / 2;
            K += Ki;
            for (j = 0; j < NBody; j++)
                if (i != j) {
                    const Body *pj = &planets[j];
                    V3 r = (pi->x - pj->x);
                    FLOAT const Uij = (pj->Gm * pi->Gm) / (r.abs());
                    U -= Uij;
                }
        }
        // remove the extra G
        return(K + U) / G;
    }



    void Predict(FLOAT dt) {
        const FLOAT dt2 = (1. / 2.) * dt;
        const FLOAT dt3 = (1. / 3.) * dt;
        int i;
        for (i = 0; i < NBody; i++) {

        }
       // xp = x + dt * (p.vel + dt2 * (p.acc + dt3 * p.jrk));
       // vel = p.vel + dt * (p.acc + dt2 * (p.jrk));
    }



    void integrate(FLOAT dt) {
        int i, j;

        // calculate predicted position for all bodies
        for (i = 0; i < NBody; i++) {
            Body *pi = &planets[i];

            V3 ai = V3(0, 0, 0);
            for (j = 0; j < NBody; j++)
                if (i != j) {
                    Body *pj = &planets[j];
                    ai += Body::acceleration(pj->x - pi->x, pj->v - pi->v, pi->Gm, pj->Gm);
                }
            pi->a = ai;
            
            //integration step:
            pi->xp = pi->x + pi->v * dt; // +pi->a * dt * dt / 2;
            pi->vp = pi->v + pi->a * dt;
            pi->tp = pi->t + dt;
        }

        // TODO(costel): calculate correction
        for (i = 0; i < NBody; i++) {
            Body *pi = &planets[i];
            pi->xp = pi->xp;
            pi->vp = pi->vp;
        }

        // Update position
        for (i = 0; i < NBody; i++) {
            Body *pi = &planets[i];
            pi->x = pi->xp;
            pi->v = pi->vp;
            pi->t = pi->tp;
        }
    }



};



typedef enum { Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto } CelestialBody;

class System {
public:
    static const int NBody = 10;
    V3 x[NBody], v[NBody], a[NBody], jerk[NBody];
    V3 xp[NBody], vp[NBody], ap[NBody], jp[NBody];
    double m[NBody];

    void init() {
        //init with JPL data
        //2440400.50 JD0 Julian date of initial state vector
        const FLOAT T0 = 2440400.50;
        CelestialBody i;
        
        i = Sun;
        x[i] = V3(4.5144118714356666407e-003, 7.2282841152065867346e-004, 2.4659100492567986271e-004);
        v[i] = V3(-2.8369446340813151639e-007, 5.1811944086463255444e-006, 2.2306588340621263489e-006);
        m[i] = 2.95912208285591102582E-4;

        i = Mercury;
        x[i] = V3(3.6030663368975339466e-001, -9.4812876741771684223e-002, -8.7466840117233140436e-002);
        v[i] = V3(3.7085069798210382186e-003, 2.4854958767430945324e-002, 1.2929109014677844626e-002);
        m[i] = 4.91254745145081175785E-11;

        i = Venus;
        x[i] = V3(6.0786466491731583464e-001, -3.5518362463675619232e-001, -1.9824142909855515515e-001);
        v[i] = V3(1.1156645711264016669e-002, 1.5494075513638794325e-002, 6.2773904546696609267e-003);
        m[i] = 7.24345620963276523095E-10;

        i = Earth;
        x[i] = V3(1.0820747754938311664e-001, -9.2711110739430602933e-001, -4.0209347855944090112e-001);
        v[i] = V3(1.6833251020051496668e-002, 1.5602036176919105255e-003, 6.7646174015847273137e-004);
        m[i] = 8.88769273403302327042E-10;

        i = Mars;
        x[i] = V3(-1.2796408611369531836e-001, -1.3262618005333617013e+000, -6.0530808652523961512e-001);
        v[i] = V3(1.4481919298277924969e-002, 8.0528538390447499843e-005, -3.5188931029397090065e-004);
        m[i] = 9.54952894222405763492E-11;

        i = Jupiter;
        x[i] = V3(-5.3896824544609061333e+000, -7.7026549518616593034e-001, -1.9866431165907522014e-001);
        v[i] = V3(1.0053452569924098185e-003, -6.5298425191689416643e-003, -2.8258787532429609536e-003);
        m[i] = 2.82534210344592625472E-7;

        i = Saturn;
        x[i] = V3(7.9527768530257360864e+000, 4.5078822184006553686e+000, 1.5201955253183338898e+000);
        v[i] = V3(-3.1594662504012930114e-003, 4.3714634278372622354e-003, 1.9441395169137103763e-003);
        m[i] = 8.45946850483065929285E-8;

        i = Uranus;
        x[i] = V3(-1.8278236586353147533e+001, -9.5764572881482056433e-001, -1.6132190397271035415e-001);
        v[i] = V3(1.7108310564806817248e-004, -3.7646704682815900043e-003, -1.6519678610257000136e-003);
        m[i] = 1.28881623813803488851E-8;

        i = Neptune;
        x[i] = V3(-1.6367191358770888335e+001, -2.3760896725373076342e+001, -9.3213866179497290101e+000);
        v[i] = V3(2.6225242764289213785e-003, -1.5277473123858904045e-003, -6.9183197562182804864e-004);
        m[i] = 1.53211248128427618918E-8;

        i = Pluto;
        x[i] = V3(-3.0447680255169362534e+001, -5.3177934960261367037e-001, 9.0596584886274922101e+000);
        v[i] = V3(2.8177758090360373050e-004, -3.1469590804946202045e-003, -1.0794238049289112837e-003);
        m[i] = 2.27624775186369921644E-12;

    }

    // given x, v calculate a, jerk
    void updateForces(V3* x, V3* v, V3* a, V3* jerk) {
        int i, j;
        V3 dx, dv, da, dj;
        double r2, ir2, ir3;
        for (i = 0; i < NBody; i++) {
            a[i] = V3(0, 0, 0);
            jerk[i] = V3(0, 0, 0);
        }

        for (i = 0; i < NBody; i++)
            for (j = i+1; j < NBody; j++) {
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
    }

    // given x, v calculate a, jerk
    void updateForcesFast(V3* x, V3* v, V3* a, V3* jerk) {
        int i, j;
        double dx, dy, dz, dvx, dvy, dvz, dax, day, daz, djx, djy, djz;
        double r2, ir2, ir3;

        for (i = 0; i < NBody; i++) {
            a[i] = V3(0, 0, 0);
            jerk[i] = V3(0, 0, 0);
        }

        for (i = 0; i < NBody; i++)
            for (j = i + 1; j < NBody; j++) {
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
    }



    void integrate(double dt) {
        int i, j;
        double dt2 = dt * dt;

        updateForcesFast(x, v, a, jerk);

        //prediction
        for (i = 0; i < NBody; i++) {
            xp[i] = x[i] + v[i] * dt + a[i] * (dt2 / 2) + jerk[i] * (dt2 * dt / 6);
            vp[i] = v[i] + a[i] * dt + jerk[i] * (dt2 / 2);
        }
        updateForcesFast(xp, vp, ap, jp);

        // correction
        for (i = 0; i < NBody; i++) {
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
        for (i = 0; i < NBody; i++) {
            double const Ki = (v[i] * v[i]) * (m[i] / 2);
            K += Ki;
            for (j = i+1; j < NBody; j++) {
                V3 r = x[i] - x[j];
                double const Uij = (m[i] * m[j]) / sqrt(r * r);
                U -= Uij;
            }
        }
        // remove the extra G (each mass is Gm)
        return(K + U)/G;
    }

    static void test() {
        System s;
        s.init();

        double dt;
        double t;

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


        double e0 = s.energy();
        printf("E(t0)=%g\n", e0);
        dt = 0.1;
        int nsteps = 0;
        __int64 t0 = now();
        for (t = 2440400.50; t <= 2440800.50; t += dt) {
            nsteps++;
            s.integrate(dt);
        }
        __int64 t1 = now();
        printf("at %g after %d steps, relative delta Energy=%g\n", 
            t, nsteps, (s.energy()-e0)/e0);

        printf("ticks for loop: %d\n", (int)((t1 - t0) / nsteps));

        { // check positional error
            V3 xtest, vtest, dx, dv;
            
            // big error for Mercury, need general relativity
            xtest = V3(3.6030663368975339466e-001, -9.4812876741771684223e-002, -8.7466840117233140436e-002);
            vtest = V3(3.7085069798210382186e-003, 2.4854958767430945324e-002, 1.2929109014677844626e-002);
            dx = s.x[Mercury] - xtest;
            printf("mercury dx=(%g %g %g)\n", dx.x, dx.y, dx.z);
            dv = s.v[Mercury] - vtest;
            printf("mercury dv=(%g %g %g)\n", dv.x, dv.y, dv.z);

            // very small error for Jupiter, newtonian aproximation works well
            xtest = V3(-4.2373615722355993645e+000, -3.1464733091277224306e+000, -1.2462241659031373286e+000);
            vtest = V3(4.6228653666202124358e-003, -5.0555918692615562352e-003, -2.2818455722493645242e-003);
            dx = s.x[Jupiter] - xtest;
            printf("jupiter dx=(%g %g %g)\n", dx.x, dx.y, dx.z);
            dv = s.v[Jupiter] - vtest;
            printf("Jupiter dv=(%g %g %g)\n", dv.x, dv.y, dv.z);

        }
    }

};


class Integrators {
public:
    V3 pos, vel;
    double mass;
    double e0;

    typedef void(Integrators::*integrate_func)(double dt);

    void init(double mass_, V3 pos_, V3 vel_) {
        mass = mass_;
        pos = pos_;
        vel = vel_;
        InitEnergy();
    };

    NOINLINE void evolve(double dt, double dt_dia, double dt_out, double dt_end, integrate_func integrate) {
        double time = 0;
        int nsteps = 0;
        double t_dia = dt_dia - dt / 2;
        double t_out = dt_out - dt / 2;
        double t_end = dt_end - dt / 2;

        __int64 t0 = now();
        while (time < t_end) {
            (this->*integrate)(dt);
            time += dt;
            nsteps++;
        }
        __int64 t1 = now();
        printf("ticks per loop:%ld\n", (int)((t1 - t0) / nsteps));
        Diagnostics(nsteps, time, dt);
    };

    V3 acc() {
        double r2 = pos * pos;
        double r3 = r2 * sqrt(r2);
        return pos * (-mass / r3);
    }

    V3 jerk() {
        double r2 = pos * pos;
        double r3 = r2 * sqrt(r2);
        return (vel + pos * (-3 * (pos * vel) / r2)) * (-mass / r3);
    }

    void updateForce(const V3& x, const V3& v, V3& a, V3& j) {
        double r2 = x * x;
        double r3 = r2 * sqrt(r2);
        a = x * (-mass / r3);
        j = (v + x * (-3 * (x * v) / r2)) * (-mass / r3);
    }


    void hermite(double dt) {
        V3 old_pos = pos;
        V3 old_vel = vel;
        V3 old_acc = acc();
        V3 old_jerk = jerk();
        double dt2 = dt * dt;
        pos += vel * dt + old_acc * (dt2 / 2) + old_jerk * (dt2 * dt / 6);
        vel += old_acc*dt + old_jerk * (dt2 / 2);
        V3 new_acc = acc();
        V3 new_jerk = jerk();
        vel = old_vel + (old_acc + new_acc) * (dt / 2) + (old_jerk - new_jerk) * (dt2 / 12);
        pos = old_pos + (old_vel + vel) * (dt / 2) + (old_acc - new_acc) * (dt2 / 12);
    }

    void hermiteRewrite(double dt) {
        V3 a, j;
        V3 ap, jp;
        V3 xp, vp;
        updateForce(pos, vel, a, j);

        //prediction
        double dt2 = dt * dt;
        xp = pos + vel * dt + a * (dt2 / 2) + j * (dt2 * dt / 6);
        vp = vel + a*dt + j * (dt2 / 2);        

        updateForce(xp, vp, ap, jp);
       
        // correction
        V3 oldv = vel;
        vel = vel + (a + ap) * (dt / 2) + (j - jp) * (dt2 / 12);
        pos = pos + (vel + oldv) * (dt / 2) + (a - ap) * (dt2 / 12);
    }


    // 90 ticks. this is an optimized version of leapfrog()
    void leapfrog2(double dt) {
        double r2, r3, f, ax, ay, az, dt2;
        double x = pos.x;
        double y = pos.y;
        double z = pos.z;
        double vx = vel.x;
        double vy = vel.y;
        double vz = vel.z;
        dt2 = 0.5 * dt;

        //acc
        r2 = x * x + y * y + z * z;
        r3 = r2 * sqrt(r2);
        f = - mass / r3;
        ax = x * f;
        ay = y * f;
        az = z * f;

        vx += ax * dt2;
        vy += ay * dt2;
        vz += az * dt2;
        
        x += vx * dt;
        y += vy * dt;
        z += vz * dt;

        //acc
        r2 = x * x + y * y + z * z;
        r3 = r2 * sqrt(r2);
        f = - mass / r3;
        ax = x * f;
        ay = y * f;
        az = z * f;

        vx += ax * dt2;
        vy += ay * dt2;
        vz += az * dt2;

        pos.x = x;
        pos.y = y;
        pos.z = z;
        vel.x = vx;
        vel.y = vy;
        vel.z = vz;
    }

    void forward(double dt) {
        V3 old_acc = acc();
        pos += vel * dt;
        vel += old_acc * dt;
    }

    // 160 ticks
    void leapfrog(double dt) {
        vel += acc() * ( dt * 0.5);
        pos += vel * dt;
        vel += acc() * (dt * 0.5);
    }

    void leapfrog1(double dt) {
        vel += acc() * (dt * 0.5);
        pos += vel * dt;
        vel += acc() * (dt * 0.5);
    }

    void hermitefast(double dt) {

        double r2, r3, f, ax, ay, az, jx, jy, jz, dt2, dt3, rv2;
        double oldx, oldy, oldz, oldvx, oldvy, oldvz;
        double oldax, olday, oldaz, oldjx, oldjy, oldjz;

        double x = pos.x;
        double y = pos.y;
        double z = pos.z;
        double vx = vel.x;
        double vy = vel.y;
        double vz = vel.z;

        //acc
        r2 = x * x + y * y + z * z;
        r3 = r2 * sqrt(r2);
        f = -mass / r3;
        ax = x * f;
        ay = y * f;
        az = z * f;

        // jer (vel + pos * (-3 * (pos * vel)/r2)) * (-mass / r3);
        rv2 = x * vx + y * vy + z * vz;
        jx = (vx + x * (-3 * rv2 / r2)) * f;
        jy = (vy + y * (-3 * rv2 / r2)) * f;
        jz = (vz + z * (-3 * rv2 / r2)) * f;

        oldx = x;
        oldy = y;
        oldz = z;
        oldvx = vx;
        oldvy = vy;
        oldvz = vz;
        oldax = ax;
        olday = ay;
        oldaz = az;
        oldjx = jx;
        oldjy = jy;
        oldjz = jz;

        dt2 = dt * dt;
        dt3 = dt2 * dt;

        x += vx * dt + ax * (dt2 / 2) + jx * (dt3 / 6);
        y += vy * dt + ay * (dt2 / 2) + jy * (dt3 / 6);
        z += vz * dt + az * (dt2 / 2) + jz * (dt3 / 6);

        vx += ax * dt + jx * (dt2 / 2);
        vy += ay * dt + jy * (dt2 / 2);
        vz += az * dt + jz * (dt2 / 2);

        //acc
        r2 = x * x + y * y + z * z;
        r3 = r2 * sqrt(r2);
        f = -mass / r3;
        ax = x * f;
        ay = y * f;
        az = z * f;

        // jerk (vel + pos * (-3 * (pos * vel)/r2)) * (-mass / r3);
        rv2 = x * vx + y * vy + z * vz;
        jx = (vx + x * (-3 * rv2 / r2)) * f;
        jy = (vy + y * (-3 * rv2 / r2)) * f;
        jz = (vz + z * (-3 * rv2 / r2)) * f;

        vx = oldvx + (oldax + ax) * (dt / 2) + (oldjx - jx) * (dt2 / 12);
        vy = oldvy + (olday + ay) * (dt / 2) + (oldjy - jy) * (dt2 / 12);
        vz = oldvz + (oldaz + az) * (dt / 2) + (oldjz - jz) * (dt2 / 12);

        x = oldx + (oldvx + vx) *(dt / 2) + (oldax - ax) * (dt2 / 12);
        y = oldy + (oldvy + vy) *(dt / 2) + (olday - ay) * (dt2 / 12);
        z = oldz + (oldvz + vz) *(dt / 2) + (oldaz - az) * (dt2 / 12);

        pos.x = x;
        pos.y = y;
        pos.z = z;
        vel.x = vx;
        vel.y = vy;
        vel.z = vz;
    }



    void rk2(double dt) {
        V3 old_pos = pos;
        V3 h_vel = vel + acc() * (0.5 * dt);
        pos += vel * (0.5 * dt);
        vel += acc() * dt;
        pos = old_pos + h_vel * dt;
    }

    void rk4(double dt) {
        double dt2 = dt * dt;
        V3 old_pos = pos;
        V3 a0 = acc();
        pos = old_pos + vel * (dt * 0.5) + a0 * 0.125 * dt2;
        V3 a1 = acc();
        pos = old_pos + vel * dt + a1 * 0.5 * dt2;
        V3 a2 = acc();
        pos = old_pos + vel * dt + (a0 + a1 * 2) * (1/6.0) * dt2;
        vel = vel + (a0 + a1 * 4 + a2) * (1 / 6.0) * dt;
    }

    double KineticEnergy() {
        return vel*vel  * (1.0 / 2);
    }

    double PotentialEnergy() {
        return -mass / sqrt(pos * pos);
    }

    void InitEnergy() {
        e0 = KineticEnergy() + PotentialEnergy();
    }

    void Diagnostics(int nsteps, double time, double dt) {
        double etot = KineticEnergy() + PotentialEnergy();
        printf("at time %g after %d steps dt=%g: \nEtot=%g delta_rel=%g\n",
            time, nsteps, dt, etot, (etot - e0) / e0);
    }

    static void test() {
        Integrators b;
        b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
        printf("e0 = %g\n", b.e0);
        double min_dt = 0.000001;
        b.evolve(1, 10, 12, 3, &Integrators::hermiteRewrite);


#if 0
        printf("Forward:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::forward);
        }

        printf("Leapfrog:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::leapfrog);
        }

        printf("Leapfrog2:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::leapfrog2);
        }

        printf("Hermite:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::hermite);
        }


#endif

        printf("Hermite rewrite:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::hermiteRewrite);
        }

#if 0

        printf("hermitefast:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::hermitefast);
        }

        printf("RK2:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::rk2);
        }

        printf("RK4:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, 100, 10, &Integrators::rk4);
        }
#endif

    }

};





int _tmain(int argc, _TCHAR* argv[])
{
#if 0
    SolarSystem ss;
    ss.init();

    printf("E(t0)=%g\n", ss.energy());

    FLOAT t;
    FLOAT dt = 0.01;
    for (t = 2440400.50; t <= 2440800.50; t += dt) {
        ss.integrate(dt);
    }
    printf("E(t1)=%g\n", ss.energy()); // energy should not change


    Body jupiter;
    jupiter.x = V3(-4.2373615722355993645e+000, -3.1464733091277224306e+000, -1.2462241659031373286e+000);
    jupiter.v = V3(4.6228653666202124358e-003, -5.0555918692615562352e-003, -2.2818455722493645242e-003);
    jupiter.t = 2440800.50;
    printf("jupiter dx=(%g %g %g)\n", jupiter.x.x - ss.planets[5].x.x, jupiter.x.y - ss.planets[5].x.y, jupiter.x.z - ss.planets[5].x.z);


    Body mercury;
    mercury.x = V3(-3.4350297408563828055e-001, -2.5118910315168408963e-001, -9.9170891227838064646e-002);
    mercury.v = V3(1.1640059959070938016e-002, -1.7964875644231269124e-002, -1.0817311005210173814e-002);
    mercury.t = 2440800.50;
    printf("mercury dx=(%g %g %g)\n", mercury.x.x - ss.planets[1].x.x, mercury.x.y - ss.planets[1].x.y, mercury.x.z - ss.planets[1].x.z);

    scanf_s("%g", t);
#else

    System::test();


//    Integrators::test();
#endif

    printf("Press a key to close.");
    getchar();

	return 0;
}

