// solar_system_simulation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include  <assert.h>
#include <cmath>
#include "vector3.h"


#define FLOAT double
#define V3 dvec3

#define G 1.0

struct Body {
    V3 x, v, a;
    FLOAT Gm;
    FLOAT t;

    // Predicted values
    V3 xp, vp;
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
            pi->xp = pi->x + pi->v * dt + pi->a * dt * dt / 2;
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




int _tmain(int argc, _TCHAR* argv[])
{
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

	return 0;
}

