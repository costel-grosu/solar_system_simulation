
/*
Tests with 1D basic integrators.
Replica of http://www.artcompsci.org/kali/development.html chapters 1,2,3.

*/

#if !defined( __BasicIntegrators_h__ )
#define __BasicIntegrators_h__

#include "Platform.h"
#include "vector3.h"

class Integrators {
public:
    V3 pos, vel;
    double mass;
    double e0;

    // integration method
    typedef void(Integrators::*integrate_func)(double dt);

    void init(double mass_, V3 pos_, V3 vel_) {
        mass = mass_;
        pos = pos_;
        vel = vel_;
        InitEnergy();
    };

    void evolve(double dt, double dt_end, integrate_func integrate) {
        double time = 0;
        int nsteps = 0;
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
        f = -mass / r3;
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
        f = -mass / r3;
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
        vel += acc() * (dt * 0.5);
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
        pos = old_pos + vel * dt + (a0 + a1 * 2) * (1 / 6.0) * dt2;
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
        b.evolve(1, 3, &Integrators::hermiteRewrite);


#if 0
        printf("Forward:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::forward);
        }

        printf("Leapfrog:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::leapfrog);
        }

        printf("Leapfrog2:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::leapfrog2);
        }

        printf("Hermite:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::hermite);
        }


#endif

        printf("Hermite rewrite:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::hermiteRewrite);
        }

#if 0

        printf("hermitefast:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::hermitefast);
        }

        printf("RK2:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::rk2);
        }

        printf("RK4:\n");
        for (double dt = 1; dt >= min_dt; dt *= 0.1) {
            b.init(1, { 1, 0, 0 }, { 0, 0.1, 0.2 });
            b.evolve(dt, 10, &Integrators::rk4);
        }
#endif

    }

};


#endif // __BasicIntegrators_h__