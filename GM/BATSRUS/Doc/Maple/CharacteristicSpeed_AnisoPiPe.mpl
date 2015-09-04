# Xing Meng
# Derivation of the wave speeds for anisotropic ion pressure + isotropic electron pressure
# 
with(Student:-VectorCalculus);
# Loading Student:-VectorCalculus
with(Student:-VectorCalculus);
with(Student:-VectorCalculus);
# 
with(Student:-VectorCalculus);
SetCoordinates('cartesian'[x, y, z]);
# 
# 
# Continuity equation:
eq1 := diff(rho(x, t), t)+Divergence(rho(x, t)*VectorField(`<,>`(ux(x, t), uy(x, t), uz(x, t)))) = 0;
# Momentum equation in X direction:
eqori2 := rho(x, t)*(diff(ux(x, t), t))+rho(x, t)*ux(x, t)*(diff(ux(x, t), x))+diff(pperp(x, t)+pe(x, t)+B(x, t)^2/(2*mu), x)+diff((Bx(x, t)*Bx(x, t))*((ppar(x, t)-pperp(x, t))/B(x, t)^2-1/mu), x) = 0;
# Expand B and divide equation by rho:
eqoori2 := subs(B(x, t) = sqrt(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), eqori2/rho(x, t));
# Eliminate dBx/dx and put back B:
eq2 := subs(diff(Bx(x, t), x) = 0, Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2 = B(x, t)^2, eqoori2);
# Momentum equation in Y direction:
eqori3 := rho(x, t)*(diff(uy(x, t), t))+rho(x, t)*ux(x, t)*(diff(uy(x, t), x))+diff(Bx(x, t)*By(x, t)*((ppar(x, t)-pperp(x, t))/B(x, t)^2-1/mu), x) = 0;
eqoori3 := subs(B(x, t) = sqrt(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), eqori3/rho(x, t));
eq3 := subs(diff(Bx(x, t), x) = 0, Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2 = B(x, t)^2, eqoori3);
# Momentum equation in Z direction:
eqori4 := rho(x, t)*(diff(uz(x, t), t))+rho(x, t)*ux(x, t)*(diff(uz(x, t), x))+diff(Bx(x, t)*Bz(x, t)*((ppar(x, t)-pperp(x, t))/B(x, t)^2-1/mu), x) = 0;
eqoori4 := subs(B(x, t) = sqrt(Bx(x, t)^2+By(x, t)^2+Bz(x, t)^2), eqori4/rho(x, t));
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
unwith(MathematicalFunctions);
# Loading LinearAlgebra


with(LinearAlgebra);
# Set ux=0 (Galilean transformation), subtract lamda*I from matrix:
CHMatrix := subs(ux(x, t) = 0, MatrixAdd(Matrix(DPMatrix), -lambda*IdentityMatrix(9)));
# Solve determinant = 0 equation for original matrix to get eigenvalues in terms of ux=-lambda, express By from B, replace Bx with cos(theta(B)).
# 
with(LinearAlgebra);
eqDP := simplify(subs(By(x, t) = sqrt(B(x, t)^2-Bx(x, t)^2-Bz(x, t)^2), Bx(x, t) = cos(theta)*B(x, t), Bx(x, t) = Bx, Bz(x, t) = Bz, B(x, t) = B, rho(x, t) = rho, ux(x, t) = ux, ppar(x, t) = ppar, pperp(x, t) = pperp, pe(x, t) = pe, pe = 3*gamma*pe*(1/5), LinearAlgebra[Determinant](Matrix(DPMatrix)) = 0));
# Divide dispersion equation by Alfven wave solution, and three entropy waves (for ppar, pperp and pe)  as a check:
eqDPs := simplify(eqDP/(ux^3*(ux^2-(B^2/mu+pperp-ppar)*cos(theta)^2/rho)));
# Fourth order equation for magnetosonic waves:
# 
collect(eqDPs, ux(x, t)^2);
# Solve the dispersion relation for all waves together:
Solution := solve(eqDP, ux);
simplify(Solution[6]);

# Loading VectorCalculus
with(VectorCalculus);
# solve for Characteristic Matrix with ux=0
eqDP2 := factor(simplify(subs(By(x, t) = sqrt(B(x, t)^2-Bx(x, t)^2-Bz(x, t)^2), Bx(x, t) = cos(theta)*B(x, t), Bx(x, t) = Bx, Bz(x, t) = Bz, B(x, t) = B, rho(x, t) = rho, ux(x, t) = ux, ppar(x, t) = ppar, pperp(x, t) = pperp, pe(x, t) = pe, pe = 3*gamma*pe*(1/5), LinearAlgebra[Determinant](Matrix(CHMatrix)) = 0)));
Solution2 := solve(eqDP2, lambda);

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
# 
# 
