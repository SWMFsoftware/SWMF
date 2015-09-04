# Xing Meng
# Derivation of the wave speeds for anisotropic ion pressure,  Semirelativistic MHD
# 
# Loading Student:-VectorCalculus
with(Student:-VectorCalculus);
with(Student:-VectorCalculus);
with(Student:-VectorCalculus);
SetCoordinates('cartesian'[x, y, z]);
# Continuity equation:
eq1 := diff(rho(x, t), t)+Divergence(rho(x, t)*VectorField(`<,>`(ux(x, t), uy(x, t), uz(x, t)))) = 0;
# Momentum equation in X direction:
eqori2 := diff(ux(x, t), t)+gamma[A]^2*`<,>`(1+V[A]^2*Bx(x, t)^2/(c^2*B(x, t)^2), V[A]^2*Bx(x, t)*By(x, t)/(c^2*B(x, t)^2), V[A]^2*Bx(x, t)*Bz(x, t)/(c^2*B(x, t)^2)).`<,>`((diff(pperp(x, t)+pe(x, t), x))/rho(x, t)+ux(x, t)*(diff(ux(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)^2/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t), ux(x, t)*(diff(uy(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)*By(x, t)/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t), ux(x, t)*(diff(uz(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)*Bz(x, t)/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t))+gamma[A]^2*`&x`(`<,>`(Bx(x, t), By(x, t), Bz(x, t)), `<,>`(0, -(diff(Bz(x, t), x)), diff(By(x, t), x))-`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(0, -(diff(Ez(x, t), x)), diff(Ey(x, t), x)))/c^2)/(mu*rho(x, t)).`<,>`(1, 0, 0) = 0;
# plug E in:

eqoori2 := simplify(subs(Ex(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(1, 0, 0), Ey(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(0, 1, 0), Ez(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(0, 0, 1), eqori2));
# Eliminate dBx/dx :
eq2 := simplify(subs(diff(Bx(x, t), x) = 0, Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2 = B(x, t)^2, eqoori2));
# Momentum equation in Y direction:
eqori3 := diff(uy(x, t), t)+gamma[A]^2*`<,>`(V[A]^2*Bx(x, t)*By(x, t)/(c^2*B(x, t)^2), 1+V[A]^2*By(x, t)/(c^2*B(x, t)^2)*By(x, t), V[A]^2*By(x, t)*Bz(x, t)/(c^2*B(x, t)^2)).`<,>`((diff(pperp(x, t)+pe(x, t), x))/rho(x, t)+ux(x, t)*(diff(ux(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)^2/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t), ux(x, t)*(diff(uy(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)*By(x, t)/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t), ux(x, t)*(diff(uz(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)*Bz(x, t)/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t))+gamma[A]^2*`&x`(`<,>`(Bx(x, t), By(x, t), Bz(x, t)), `<,>`(0, -(diff(Bz(x, t), x)), diff(By(x, t), x))-`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(0, -(diff(Ez(x, t), x)), diff(Ey(x, t), x)))/c^2)/(mu*rho(x, t)).`<,>`(0, 1, 0) = 0;
eqoori3 := simplify(subs(Ex(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(1, 0, 0), Ey(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(0, 1, 0), Ez(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(0, 0, 1), eqori3));
eq3 := subs(diff(Bx(x, t), x) = 0, eqoori3);
# Momentum equation in Z direction:
eqori4 := diff(uz(x, t), t)+gamma[A]^2*`<,>`(V[A]^2*Bx(x, t)*Bz(x, t)/(c^2*B(x, t)^2), V[A]^2*Bz(x, t)*By(x, t)/(c^2*B(x, t)^2), 1+V[A]^2*Bz(x, t)/(c^2*B(x, t)^2)*Bz(x, t)).`<,>`((diff(pperp(x, t)+pe(x, t), x))/rho(x, t)+ux(x, t)*(diff(ux(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)^2/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t), ux(x, t)*(diff(uy(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)*By(x, t)/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t), ux(x, t)*(diff(uz(x, t), x))+(diff((ppar(x, t)-pperp(x, t))*Bx(x, t)*Bz(x, t)/(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), x))/rho(x, t))+gamma[A]^2*`&x`(`<,>`(Bx(x, t), By(x, t), Bz(x, t)), `<,>`(0, -(diff(Bz(x, t), x)), diff(By(x, t), x))-`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(0, -(diff(Ez(x, t), x)), diff(Ey(x, t), x)))/c^2)/(mu*rho(x, t)).`<,>`(0, 0, 1) = 0;

eqoori4 := simplify(subs(E(x, t) = sqrt(Ex(x, t)^2+Ey(x, t)^2+Ez(x, t)^2), Ex(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(1, 0, 0), Ey(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(0, 1, 0), Ez(x, t) = -`&x`(`<,>`(ux(x, t), uy(x, t), uz(x, t)), `<,>`(Bx(x, t), By(x, t), Bz(x, t))).`<,>`(0, 0, 1), eqori4));
eq4 := subs(diff(Bx(x, t), x) = 0, Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2 = B(x, t)^2, eqoori4);

# Check that grad Pe is recognized:

Gradient(pe(x, t));
# Induction equation in X direction:


eq5 := diff(Bx(x, t), t)+diff(uex(x, t)*Bx(x, t)-uex(x, t)*Bx(x, t), x) = 0;

# Induction equation in Y direction:
eqori6 := diff(By(x, t), t)+diff(uex(x, t)*By(x, t)-Bx(x, t)*uey(x, t), x) = 0;
# Set electron velocity as uey = uy - Jy/(ne), use alpha = e/Mi (but this would result in whistler waves!):
eq6TBD := subs(diff(Bx(x, t), x) = 0, uex(x, t) = ux(x, t), uey(x, t) = uy(x, t)+(diff(Bz(x, t), x))/(alpha*mu*rho(x, t)), eqori6);
# Set electron velocity = ion velocity in the induction equation:
eq6 := subs(diff(Bx(x, t), x) = 0, uex(x, t) = ux(x, t), uey(x, t) = uy(x, t), eqori6);
eqori7 := diff(Bz(x, t), t)+diff(uex(x, t)*Bz(x, t)-Bx(x, t)*uez(x, t), x) = 0;
eq7TBD := subs(diff(Bx(x, t), x) = 0, uex(x, t) = ux(x, t), uez(x, t) = uz(x, t)-(diff(By(x, t), x))/(alpha*mu*rho(x, t)), eqori7);
eq7 := subs(diff(Bx(x, t), x) = 0, uex(x, t) = ux(x, t), uez(x, t) = uz(x, t), eqori7);
# Parallel ion pressure equation:
eq8 := diff(ppar(x, t), t)+diff(ppar(x, t)*ux(x, t), x)+2*ppar(x, t)*Bx(x, t)*(Bx(x, t)*(diff(ux(x, t), x))+By(x, t)*(diff(uy(x, t), x))+Bz(x, t)*(diff(uz(x, t), x)))/B(x, t)^2 = 0;
# Perpendicular ion pressure equation:
eq9 := diff(pperp(x, t), t)+diff(pperp(x, t)*ux(x, t), x)+pperp(x, t)*(diff(ux(x, t), x))-pperp(x, t)*Bx(x, t)*(Bx(x, t)*(diff(ux(x, t), x))+By(x, t)*(diff(uy(x, t), x))+Bz(x, t)*(diff(uz(x, t), x)))/B(x, t)^2 = 0;
# Electron pressure equation:
eqori10 := diff(pe(x, t), t)+diff(pe(x, t)*uex(x, t), x)+(2/3)*pe(x, t)*(diff(uex(x, t), x)) = 0;
eq10 := subs(uex(x, t) = ux(x, t), eqori10);
# Check VectorField() and Curl() functions:
VectorField(`<,>`(uex(x, t), uey(x, t), uez(x, t))) = VectorField(`<,>`(ux(x, t), uy(x, t), uz(x, t)))-(Curl(VectorField(`<,>`(Bx(x, t), By(x, t), Bz(x, t)))))/(alpha*mu*rho(x, t));
# Dispersion matrix = coefficients of d(rho)/dx, ... , dPe/dx:

with(Student:-VectorCalculus);
DPMatrix := linalg[genmatrix]([eq1, eq2, eq3, eq4, eq6, eq7, eq8, eq9, eq10], [diff(rho(x, t), x), diff(ux(x, t), x), diff(uy(x, t), x), diff(uz(x, t), x), diff(By(x, t), x), diff(Bz(x, t), x), diff(ppar(x, t), x), diff(pperp(x, t), x), diff(pe(x, t), x)]);
# 
# Loading LinearAlgebra


with(LinearAlgebra);
# Loading LinearAlgebra
DPMatrix00 := simplify(subs(rho(x, t) = rho, ux(x, t) = ux, uy(x, t) = uy, uz(x, t) = uz, Bx(x, t) = Bx, By(x, t) = By, Bz(x, t) = Bz, B(x, t) = B, ppar(x, t) = ppar, pperp(x, t) = pperp, pe(x, t) = pe, Matrix(DPMatrix)));
with(LinearAlgebra);
# 
# 
# parallel magnetic field
with(LinearAlgebra);
DPMatrix1 := simplify(subs(B = Bx, By = 0, Bz = 0, Bx = sqrt(rho*mu)*V[A], gamma[A]^2 = 1/(1+V[A]^2/c^2), MatrixAdd(Matrix(DPMatrix00), -lambda*IdentityMatrix(9))));
eqDP1 := simplify(LinearAlgebra[Determinant](Matrix(DPMatrix1)) = 0);
solve(eqDP1, lambda);
simplify(subs(c^2+V[A]^2 = c^2/gamma[A]^2, (rho*ux*c+sqrt(rho*c^2*pperp-rho*c^2*ppar+rho^2*c^2*V[A]^2-rho^2*V[A]^2*ux^2+rho*V[A]^2*pperp-rho*V[A]^2*ppar+rho^2*V[A]^4))*c/(rho*(c^2+V[A]^2))));

# perpendicular magnetic field
with(LinearAlgebra);
DPMatrix2 := simplify(subs(Bx = 0, B = sqrt(By^2+Bz^2), By = sqrt(V[A]^2*mu*rho-Bz^2), gamma[A]^2 = 1/(1+V[A]^2/c^2), MatrixAdd(Matrix(DPMatrix00), -lambda*IdentityMatrix(9))));
eqDP2 := simplify(LinearAlgebra[Determinant](Matrix(DPMatrix2)) = 0);
soln2 := solve(eqDP2, lambda);
simplify(subs(c^2+V[A]^2 = c^2/gamma[A]^2, V[A] = sqrt(B^2/(mu*rho)), pe = 3*gamma*pe*(1/5), soln2[9]));
# general case
# 
# 
DPMatrix0 := MatrixAdd(Matrix(DPMatrix00), -lambda*IdentityMatrix(9));
eqDP0 := simplify(subs(By = sqrt(B^2-Bx^2-Bz^2), gamma[A] = 1/sqrt(1+V[A]^2/c^2), V[A] = sqrt(B^2/(mu*rho)), LinearAlgebra[Determinant](Matrix(DPMatrix0)) = 0));
eqDP0new := factor(eqDP0);
# 
# to obtain similar form of P4(lambda) as isotropic MHD in the paper [ eqn (46) ]
collect(simplify(subs(pe = 3*gamma*pe*(1/5), 1/3*(-12*lambda*Bx^2*mu*B^2*c^2*ppar*ux*rho+6*lambda*Bx^2*mu*B^2*c^2*pperp*ux*rho+6*Bx^2*ux^2*ppar*B^2*rho*mu*c^2-3*Bx^2*mu*B^2*c^2*pperp*ux^2*rho+6*lambda^2*Bx^2*mu*B^2*c^2*ppar*rho-3*lambda^2*Bx^2*mu*B^2*c^2*pperp*rho-12*lambda*B^4*rho*c^2*mu*ux*pperp-10*lambda*B^4*pe*ux*c^2*rho*mu+3*B^6*c^2*ux^2*rho+9*Bx^4*c^2*ppar^2*mu-5*B^4*Bx^2*pe*c^2-9*B^4*Bx^2*c^2*ppar-3*Bx^4*pperp^2*c^2*mu+3*lambda^2*B^6*c^2*rho-15*B^2*Bx^2*pe*ppar*c^2*mu+5*B^4*pe*ux^2*c^2*rho*mu-18*B^2*Bx^2*ppar*c^2*pperp*mu+12*lambda*ux^3*B^4*rho^2*mu*c^2+6*B^4*ux^2*c^2*pperp*rho*mu+6*lambda^2*B^4*rho*c^2*mu*pperp-18*lambda^2*ux^2*B^4*rho^2*mu*c^2+12*lambda^3*B^4*rho^2*ux*mu*c^2+5*lambda^2*B^4*pe*rho*c^2*mu-3*B^6*lambda^4*rho-5*Bx^4*pe*c^2*pperp*mu-3*ux^4*B^4*rho^2*mu*c^2+3*B^2*Bx^2*pperp^2*c^2*mu+9*Bx^4*c^2*pperp*ppar*mu+20*Bx^4*pe*c^2*ppar*mu-3*lambda^4*B^4*rho^2*c^2*mu-6*lambda*B^6*c^2*ux*rho+9*B^4*Bx^2*lambda^2*ppar+6*B^6*lambda^3*ux*rho-3*B^6*lambda^2*ux^2*rho+5*B^4*Bx^2*pe*lambda^2))+simplify(subs(V[A]^2 = B^2/(mu*rho), a^2 = (3*ppar+gamma*pe)/rho, bx = Bx/B, B^4*rho^2*c^2*mu*((lambda-ux)^4-a^2*(lambda-ux)^2-(c^2-lambda^2)*V[A]^2*((lambda-ux)^2-a^2*bx^2)/c^2)))), lambda);
# compare with AnisoP
(2*B^4*pperp*c^2*rho*mu-3*B^4*ppar*mu*rho*c^2+2*B^2*Bx^2*ppar*rho*mu*c^2-B^2*Bx^2*pperp*c^2*rho*mu)*lambda^2+(-4*B^2*Bx^2*ux*ppar*rho*mu*c^2+2*B^2*Bx^2*ux*c^2*pperp*rho*mu-4*B^4*pperp*ux*c^2*rho*mu+6*B^4*ppar*mu*rho*c^2*ux)*lambda+3*Bx^4*c^2*pperp*ppar*mu-6*B^2*Bx^2*ppar*c^2*pperp*mu+2*B^4*ux^2*c^2*pperp*rho*mu+3*Bx^4*c^2*ppar^2*mu-3*B^4*ppar*rho*mu*c^2*ux^2+2*Bx^2*ux^2*ppar*B^2*rho*mu*c^2-Bx^2*mu*B^2*c^2*pperp*ux^2*rho+B^2*Bx^2*pperp^2*c^2*mu-Bx^4*pperp^2*c^2*mu-(2*B^4*pperp*c^2*rho*mu-3*B^4*ppar*mu*rho*c^2+2*B^2*Bx^2*ppar*rho*mu*c^2-B^2*Bx^2*pperp*c^2*rho*mu)*lambda^2-(-4*B^2*Bx^2*ux*ppar*rho*mu*c^2+2*B^2*Bx^2*ux*c^2*pperp*rho*mu-4*B^4*pperp*ux*c^2*rho*mu+6*B^4*ppar*mu*rho*c^2*ux)*lambda-3*Bx^4*c^2*ppar^2*mu+Bx^4*pperp^2*c^2*mu+Bx^2*mu*B^2*c^2*pperp*ux^2*rho-2*Bx^2*ux^2*ppar*B^2*rho*mu*c^2+3*B^2*Bx^2*gamma*pe*ppar*c^2*mu-2*B^4*ux^2*c^2*pperp*rho*mu-4*Bx^4*gamma*pe*c^2*ppar*mu+6*B^2*Bx^2*ppar*c^2*pperp*mu-3*Bx^4*c^2*pperp*ppar*mu+3*B^4*ppar*rho*mu*c^2*ux^2+Bx^4*gamma*pe*c^2*pperp*mu-B^2*Bx^2*pperp^2*c^2*mu;
# see latex equation file for final form of the dispersion relation
# 
soln := solve(eqDP0, lambda);
# Alfven waves, soln[1] and soln[2]
simplify(subs(sqrt(B^2-Bx^2-Bz^2) = B*byy, Bx = B*bx, Bz = B*bz, soln[2]));
# special case: zero flow velocity
eqDP0noflow := simplify(subs(ux = 0, uy = 0, uz = 0, eqDP0new));
solnnoflow := solve(eqDP0noflow, lambda);
# 
factor(simplify(subs(B^2+c^2*rho*mu = c^2*rho*mu/gamma[A]^2, Bx = cos(theta)*B, pe = 3*gamma*pe*(1/5), eqDP0noflow)));
# 
# 
# fast and slow magnetosonic waves for zero flow velocity
# 
simplify(subs(B^2+c^2*rho*mu = c^2*rho*mu/gamma[A]^2, Bx = cos(theta)*B, pe = 3*gamma*pe*(1/5), solnnoflow[4]));
# 
# 
# 

fast := simplify(subs(collect((B^4*gamma*pe*cos(theta)^2+2*cos(theta)^2*B^2*mu*c^2*ppar*rho-cos(theta)^2*B^2*mu*c^2*pperp*rho+B^4*rho*c^2+3*B^4*ppar*cos(theta)^2+B^2*gamma*pe*rho*c^2*mu+2*B^2*rho*c^2*mu*pperp)/B^2, B)));
factor(simplify(9*B^4*ppar^2*cos(theta)^4+B^4*c^4*rho^2+4*B^2*c^2*rho*mu*pperp*gamma*pe*cos(theta)^2-fast^2+2*B^2*gamma^2*pe^2*c^2*rho*mu*cos(theta)^2-2*cos(theta)^2*mu^2*c^4*pperp*rho^2*gamma*pe-12*B^2*c^2*rho*mu*pperp*ppar*cos(theta)^2-4*rho^2*c^4*mu*B^2*cos(theta)^2*gamma*pe-8*cos(theta)^2*ppar*rho^2*mu^2*c^4*gamma*pe+6*cos(theta)^4*B^2*mu*c^2*pperp*rho*ppar-4*rho^2*c^4*mu^2*cos(theta)^4*gamma*pe*pperp+16*rho^2*c^4*mu^2*cos(theta)^4*gamma*pe*ppar-6*cos(theta)^4*B^2*mu*c^2*pperp*rho*gamma*pe-6*B^2*gamma*pe*c^2*rho*mu*ppar*cos(theta)^2+20*cos(theta)^4*B^2*ppar*rho*mu*c^2*gamma*pe+6*B^4*gamma*pe*cos(theta)^4*ppar-6*B^4*ppar*cos(theta)^2*rho*c^2-3*cos(theta)^4*mu^2*c^4*pperp^2*rho^2+16*cos(theta)^4*ppar^2*rho^2*mu^2*c^4+gamma^2*pe^2*c^4*rho^2*mu^2+4*B^2*c^4*rho^2*mu*pperp+B^4*gamma^2*pe^2*cos(theta)^4+4*c^4*rho^2*mu^2*pperp^2-4*rho*B^2*cos(theta)^4*pperp^2*c^2*mu+4*rho*B^2*cos(theta)^2*pperp^2*c^2*mu-2*cos(theta)^2*B^2*mu*c^4*pperp*rho^2-16*cos(theta)^2*ppar*rho^2*mu^2*c^4*pperp+8*cos(theta)^4*mu^2*c^4*pperp*rho^2*ppar+4*c^4*rho^2*mu^2*pperp*gamma*pe-2*B^4*gamma*pe*cos(theta)^2*rho*c^2+24*cos(theta)^4*B^2*ppar^2*rho*mu*c^2-8*cos(theta)^2*B^2*ppar*rho^2*mu*c^4+2*B^2*gamma*pe*c^4*rho^2*mu));
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
