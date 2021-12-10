#![allow(non_snake_case)]

pub use anyhow::Error;
use anyhow::{bail, ensure, Result};
pub use ndarray::prelude::*;
pub use ndarray_linalg::cholesky::*;
use std::convert::TryInto;

// pub enum Error {}

type Matrix = Array2<f64>;
// type MatrixViewMut = ArrayViewMut2<f64>;
type Vector = Array1<f64>;
type IVector = Array1<i32>;
type BVector = Array1<bool>;

pub fn solve_quadprog(
    mut G: Matrix,
    g0: Vector,
    CE: Matrix,
    ce0: Vector,
    CI: Matrix,
    ci0: Vector,
    x: &mut Vector,
) -> Result<f64> {
    let n = G.ncols();
    let p = CE.ncols();
    let m = CI.ncols();

    ensure!(G.nrows() == n);
    ensure!(CE.nrows() == n);
    ensure!(CI.nrows() == n);
    ensure!(ce0.len() == p);
    ensure!(ci0.len() == m);
    ensure!(x.len() == n);
    // let mut x: Vector = Vector::zeros(n);

    let mut R = Matrix::zeros((n, n));

    let mut s = Vector::zeros(m + p);
    let mut z = Vector::zeros(n);
    let mut r = Vector::zeros(m + p);
    let mut d = Vector::zeros(n);
    let mut np = Vector::zeros(n);
    let mut u = Vector::zeros(m + p);
    let mut x_old = Vector::zeros(n);
    let mut u_old = Vector::zeros(m + p);
    let mut J = Matrix::zeros((n, n));

    let mut f_value: f64;
    let mut psi: f64;
    let mut sum: f64;
    let mut ss: f64;

    let inf = f64::INFINITY;
    // let inf: f64;
    // if std::numeric_limits<double>::has_infinity
    //   inf = std::numeric_limits<double>::infinity();
    // else
    //   inf = 1.0E300;

    /* t is the step lenght, which is the minimum of the partial step length t1
     * and the full step length t2 */
    let mut t: f64;
    let mut t1: f64;
    let mut t2: f64;
    let mut A = IVector::zeros(m + p);
    let mut A_old = IVector::zeros(m + p);
    let mut iai = IVector::zeros(m + p);
    let mut iq: usize;
    // TODO cap upper limit of iterations?
    // let mut iter: usize = 0;
    let mut iaexcl = BVector::default(m + p);

    /*
     * Preprocessing phase
     */

    /* compute the trace of the original matrix G */
    let c1: f64 = G.diag().sum();
    /* decompose the matrix G in the form L^T L */
    // TODO it says L in the comment, but the function uses the upper triangle.

    cholesky_decomposition(&mut G)?;
    // G.cholesky_inplace(UPLO::Upper)?;

    /* initialize the matrix R */
    let mut R_norm = 1.0; /* this variable will hold the norm of the matrix R */

    /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
    let mut c2 = 0.0;
    for i in 0..n {
        // d should be a vector with 1 in the ith position and zeros
        // elsewhere. This is just a performance oriented construction.
        d[i] = 1.0;
        forward_elimination(&G, &mut z, &d);
        for j in 0..n {
            J[(i, j)] = z[j];
        }
        c2 += z[i];
        d[i] = 0.0;
    }

    /* c1 * c2 is an estimate for cond(G) */

    /*
     * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x
     * this is a feasible point in the dual space
     * x = G^-1 * g0
     */
    cholesky_solve(&G, x, &g0);
    for i in 0..n {
        x[i] = -x[i];
    }
    /* and compute the current solution value */
    f_value = 0.5 * g0.dot(x);

    /* Add equality constraints to the working set A */
    iq = 0;
    for i in 0..(p) {
        for j in 0..(n) {
            np[j] = CE[(j, i)];
        }
        compute_d(&mut d, &J, &np);
        update_z(&mut z, &J, &d, iq);
        update_r(&R, &mut r, &d, iq);

        /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint
        becomes feasible */
        t2 = 0.0;
        // i.e. z != 0
        if z.dot(&z) > f64::EPSILON {
            t2 = (-np.dot(x) - ce0[i]) / z.dot(&np);
        }

        /* set x = x + t2 * z */
        for k in 0..(n) {
            x[k] += t2 * z[k];
        }

        /* set u = u+ */
        u[iq] = t2;
        for k in 0..(iq) {
            u[k] -= t2 * r[k];
        }

        /* compute the new solution value */
        f_value += 0.5 * (t2 * t2) * z.dot(&np);
        // TODO why is this negative?
        A[i] = -(i as i32) - 1;

        if !add_constraint(&mut R, &mut J, &mut d, &mut iq, &mut R_norm) {
            // Equality constraints are linearly dependent
            bail!("Constraints are linearly dependent");
            // return f_value;
        }
    }

    /* set iai = K \ A */
    for i in (0)..(m) {
        // TODO i32?
        iai[i] = i as i32;
    }

    let mut ip: usize; // this is the index of the constraint to be added to the active set
    'l1: loop {
        // iter += 1;
        /* step 1: choose a violated constraint */
        for i in (p)..(iq) {
            ip = A[i].try_into()?;
            iai[ip] = -1;
        }

        /* compute s[x] = ci^T * x + ci0 for all elements of K \ A */
        ss = 0.0;
        psi = 0.0; /* this value will contain the sum of all infeasibilities */
        ip = 0; /* ip will be the index of the chosen violated constraint */
        for i in (0)..(m) {
            iaexcl[i] = true;
            sum = 0.0;
            for j in (0)..(n) {
                sum += CI[(j, i)] * x[j];
            }
            sum += ci0[i];
            s[i] = sum;
            psi += 0.0_f64.min(sum);
        }

        if (psi).abs() <= m as f64 * f64::EPSILON * c1 * c2 * 100.0 {
            /* numerically there are not infeasibilities anymore */
            return Ok(f_value);
        }

        /* save old values for u and A */
        for i in (0)..(iq) {
            u_old[i] = u[i];
            A_old[i] = A[i];
        }
        /* and for x */
        for i in (0)..(n) {
            x_old[i] = x[i];
        }

        'l2: loop {
            /* Step 2: check for feasibility and determine a new S-pair */
            for i in (0)..(m) {
                if s[i] < ss && iai[i] != -1 && iaexcl[i] {
                    ss = s[i];
                    ip = i;
                }
            }
            if ss >= 0.0 {
                return Ok(f_value);
            }

            /* set np = n[ip] */
            for i in (0)..(n) {
                np[i] = CI[(i, ip)];
            }
            /* set u = [u 0]^T */
            u[iq] = 0.0;
            /* add ip to the active set A */
            A[iq] = ip as i32;

            'l2a: loop {
                /* Step 2a: determine step direction */
                /* compute z = H np: the step direction in the primal space (through J, see the paper) */
                compute_d(&mut d, &J, &np);
                update_z(&mut z, &J, &d, iq);
                /* compute N* np (if q > 0): the negative of the step direction in the dual space */
                update_r(&R, &mut r, &d, iq);

                /* Step 2b: compute step length */
                let mut l: usize = 0;
                /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
                t1 = inf; /* +inf */
                /* find the index l s.t. it reaches the minimum of u+[x] / r */
                for k in (p)..(iq) {
                    if r[k] > 0.0 {
                        if u[k] / r[k] < t1 {
                            t1 = u[k] / r[k];
                            // TODO... this can be negative.
                            l = A[k].try_into()?;
                        }
                    }
                }
                /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
                if (z.dot(&z).abs()) > f64::EPSILON
                // i.e. z != 0
                {
                    t2 = -s[ip] / z.dot(&np);
                    if t2 < 0.0_f64
                    // patch suggested by Takano Akio for handling numerical inconsistencies
                    {
                        t2 = inf;
                    }
                } else {
                    t2 = inf; /* +inf */
                }

                /* the step is chosen as the minimum of t1 and t2 */
                t = t1.min(t2);

                /* Step 2c: determine new S-pair and take step: */

                /* case (i): no step in primal or dual space */
                if t >= inf {
                    /* QPP is infeasible */
                    // FIXME: unbounded to raise
                    return Ok(inf);
                }
                /* case (ii): step in dual space */
                if t2 >= inf {
                    /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
                    for k in (0)..(iq) {
                        u[k] -= t * r[k];
                    }
                    u[iq] += t;
                    iai[l] = l as i32;
                    delete_constraint(&mut R, &mut J, &mut A, &mut u, n, p, &mut iq, l)?;
                    continue 'l2a;
                }

                /* case (iii): step in primal and dual space */

                /* set x = x + t * z */
                for k in (0)..(n) {
                    x[k] += t * z[k];
                }
                /* update the solution value */
                f_value += t * z.dot(&np) * (0.5 * t + u[iq]);
                /* u = u + t * [-r 1] */
                for k in (0)..(iq) {
                    u[k] -= t * r[k];
                }
                u[iq] += t;

                if (t - t2).abs() < f64::EPSILON {
                    /* full step has taken */
                    /* add constraint ip to the active set*/
                    if !add_constraint(&mut R, &mut J, &mut d, &mut iq, &mut R_norm) {
                        iaexcl[ip] = false;
                        delete_constraint(&mut R, &mut J, &mut A, &mut u, n, p, &mut iq, ip)?;
                        for i in (0)..(m) {
                            iai[i] = i as i32;
                        }
                        for i in (p)..(iq) {
                            A[i] = A_old[i];
                            u[i] = u_old[i];
                            // TODO.. use try_into() instead? but hot loop...
                            iai[A[i] as usize] = -1;
                        }
                        for i in (0)..(n) {
                            x[i] = x_old[i];
                        }
                        continue 'l2; /* go to step 2 */
                    } else {
                        iai[ip] = -1;
                    }
                    continue 'l1;
                }

                /* a patial step has taken */
                /* drop constraint l */
                iai[l] = l as i32;
                delete_constraint(&mut R, &mut J, &mut A, &mut u, n, p, &mut iq, l)?;

                /* update s[ip] = CI * x + ci0 */
                sum = 0.0;
                for k in (0)..(n) {
                    sum += CI[(k, ip)] * x[k];
                }
                s[ip] = sum + ci0[ip];
            }
        }
    }
}

fn distance(a: f64, b: f64) -> f64 {
    let a1 = (a).abs();
    let b1 = (b).abs();
    if a1 > b1 {
        let t = b1 / a1;
        return a1 * (1.0 + t * t).sqrt();
    } else if b1 > a1 {
        let t = a1 / b1;
        return b1 * (1.0 + t * t).sqrt();
    }
    return a1 * 2.0_f64.sqrt();
}

fn add_constraint(
    R: &mut Matrix,
    J: &mut Matrix,
    d: &mut Vector,
    iq: &mut usize,
    R_norm: &mut f64,
) -> bool {
    let n = d.len();
    let mut cc: f64;
    let mut ss: f64;
    let mut h: f64;
    let mut t1: f64;
    let mut t2: f64;
    let mut xny: f64;

    /* we have to find the Givens rotation which will reduce the element
    d[j] to zero.
    if it is already zero we don't have to do anything, except of
    decreasing j */
    for j in ((*iq + 1)..n).rev() {
        /* The Givens rotation is done with the matrix (cc cs, cs -cc).
        If cc is one, then element (j) of d is zero compared with element
        (j - 1). Hence we don't have to do anything.
        If cc is zero, then we just have to switch column (j) and column (j - 1)
        of J. Since we only switch columns in J, we have to be careful how we
        update d depending on the sign of gs.
        Otherwise we have to apply the Givens rotation to these columns.
        The i - 1 element of d has to be updated to h. */
        cc = d[j - 1];
        ss = d[j];
        h = distance(cc, ss);
        if h.abs() < f64::EPSILON
        // h == 0
        {
            continue;
        }
        d[j] = 0.0;
        ss = ss / h;
        cc = cc / h;
        if cc < 0.0 {
            cc = -cc;
            ss = -ss;
            d[j - 1] = -h;
        } else {
            d[j - 1] = h;
        }
        xny = ss / (1.0 + cc);
        for k in (0)..(n) {
            t1 = J[(k, j - 1)];
            t2 = J[(k, j)];
            J[(k, j - 1)] = t1 * cc + t2 * ss;
            J[(k, j)] = xny * (t1 + J[(k, j - 1)]) - t2;
        }
    }
    /* update the number of constraints added*/
    *iq += 1;
    /* To update R we have to put the iq components of the d vector
    into column iq - 1 of R
    */
    for i in (0)..(*iq) {
        R[(i, *iq - 1)] = d[i];
    }

    if (d[*iq - 1]).abs() <= f64::EPSILON * *R_norm {
        // problem degenerate
        return false;
    }
    *R_norm = R_norm.max((d[*iq - 1]).abs());
    return true;
}

fn forward_elimination(L: &Matrix, y: &mut Vector, b: &Vector) {
    let n = L.nrows();

    y[0] = b[0] / L[(0, 0)];
    for i in 1..n {
        y[i] = b[i];
        for j in 0..i {
            y[i] -= L[(i, j)] * y[j];
        }
        y[i] = y[i] / L[(i, i)];
    }
}

fn cholesky_solve(L: &Matrix, x: &mut Vector, b: &Vector) {
    let n = L.nrows();
    // TODO @stash
    let mut y = Vector::zeros(n);

    /* Solve L * y = b */
    forward_elimination(L, &mut y, b);
    /* Solve L^T * x = y */
    backward_elimination(L, x, &y);
}

fn backward_elimination(U: &Matrix, x: &mut Vector, y: &Vector) {
    let n = U.nrows();

    x[n - 1] = y[n - 1] / U[(n - 1, n - 1)];
    // TODO (i = n - 2; i >= 0; i--)
    for i in (0..(n - 1)).rev() {
        x[i] = y[i];
        for j in (i + 1)..n {
            x[i] -= U[(i, j)] * x[j];
        }
        x[i] = x[i] / U[(i, i)];
    }
}

fn compute_d(d: &mut Vector, J: &Matrix, np: &Vector) {
    let n = d.len();

    /* compute d = H^T * np */
    for i in 0..(n) {
        let mut sum: f64 = 0.0f64;
        for j in 0..(n) {
            sum += J[(j, i)] * np[j];
        }
        d[i] = sum;
    }
}

fn update_z(z: &mut Vector, J: &Matrix, d: &Vector, iq: usize) {
    let n = z.len();

    /* setting of z = H * d */
    for i in 0..(n) {
        z[i] = 0.0;
        for j in (iq)..(n) {
            z[i] += J[(i, j)] * d[j];
        }
    }
}

fn update_r(R: &Matrix, r: &mut Vector, d: &Vector, iq: usize) {
    /* setting of r = R^-1 d */
    for i in (0..(iq)).rev() {
        let mut sum: f64 = 0.0f64;
        for j in (i + 1)..(iq) {
            sum += R[(i, j)] * r[j];
        }
        r[i] = (d[i] - sum) / R[(i, i)];
    }
}

fn delete_constraint(
    R: &mut Matrix,
    J: &mut Matrix,
    A: &mut IVector,
    u: &mut Vector,
    n: usize,
    p: usize,
    iq: &mut usize,
    l: usize,
) -> Result<()> {
    let mut qq: usize = 0; // just to prevent warnings from smart compilers
    let mut cc: f64;
    let mut ss: f64;
    let mut h: f64;
    let mut xny: f64;
    let mut t1: f64;
    let mut t2: f64;

    let mut found = false;
    /* Find the index qq for active constraint l to be removed */
    for i in (p)..(*iq) {
        if A[i] == l as i32 {
            qq = i;
            found = true;
            break;
        }
    }

    if !found {
        bail!(
            "Attempt to delete non existing constraint, constraint: {:?}",
            l
        );
    }
    /* remove the constraint from the active set and the duals */
    for i in (qq)..(*iq - 1) {
        A[i] = A[i + 1];
        u[i] = u[i + 1];
        for j in (0)..(n) {
            R[(j, i)] = R[(j, i + 1)];
        }
    }

    A[*iq - 1] = A[*iq];
    u[*iq - 1] = u[*iq];
    A[*iq] = 0;
    u[*iq] = 0.0;
    for j in (0)..(*iq) {
        R[(j, *iq - 1)] = 0.0;
    }
    /* constraint has been fully removed */
    *iq -= 1;

    if *iq == 0 {
        return Ok(());
    }

    for j in (qq)..(*iq) {
        cc = R[(j, j)];
        ss = R[(j + 1, j)];
        h = distance(cc, ss);
        if h.abs() < f64::EPSILON
        // h == 0
        {
            continue;
        }
        cc = cc / h;
        ss = ss / h;
        R[(j + 1, j)] = 0.0;
        if cc < 0.0 {
            R[(j, j)] = -h;
            cc = -cc;
            ss = -ss;
        } else {
            R[(j, j)] = h;
        }

        xny = ss / (1.0 + cc);
        for k in (j + 1)..(*iq) {
            t1 = R[(j, k)];
            t2 = R[(j + 1, k)];
            R[(j, k)] = t1 * cc + t2 * ss;
            R[(j + 1, k)] = xny * (t1 + R[(j, k)]) - t2;
        }
        for k in (0)..(n) {
            t1 = J[(k, j)];
            t2 = J[(k, j + 1)];
            J[(k, j)] = t1 * cc + t2 * ss;
            J[(k, j + 1)] = xny * (J[(k, j)] + t1) - t2;
        }
    }
    Ok(())
}

fn cholesky_decomposition(A: &mut Matrix) -> Result<()> {
    let n = A.nrows();

    for i in (0)..(n) {
        for j in (i)..(n) {
            let mut sum: f64 = A[(i, j)];
            for k in (0..i).rev() {
                sum -= A[(i, k)] * A[(j, k)];
            }
            if i == j {
                if sum <= 0.0 {
                    bail!("Error in cholesky decomposition\nsum: {}\nA: {:?}", sum, A);
                }
                A[(i, i)] = (sum).sqrt();
            } else {
                A[(j, i)] = sum / A[(i, i)];
            }
        }
        for k in (i + 1)..(n) {
            A[(i, k)] = A[(k, i)];
        }
    }
    Ok(())
}

/*
%s/const Vector<double>& \([a-zA-Z_]\+\)/\1: \&Vector
%s/const Matrix<double>& \([a-zA-Z_]\+\)/\1: \&Matrix
%s/Vector<double>& \([a-zA-Z_]\+\)/\1: \&mut Vector
%s/\[\([^\]]\+\)\]\[\([^\]]\+\)\]/[(\1, \2)]/g
%s/Vector<double> \(\(\w\+\)(\([^)]\+\)),\?\s*;\?\)\+/let \1 = Vector::zeros(\2);/g
*/

#[cfg(test)]
mod tests {
    use approx::{assert_abs_diff_eq, assert_ulps_eq};
    use ndarray::{array, Array, Array2};

    use super::*;

    // #[test]
    // fn cholesky_decom_is_legit() -> Result<()> {
    //     #[rustfmt::skip]
    //     let g = array![
    //         [4.0, -2.0],
    //         [-2.0, 4.0],
    //     ];
    //     let mut manual = g.clone();
    //     let mut linalg = g;
    //     cholesky_decomposition(&mut manual)?;
    //     linalg.cholesky_inplace(UPLO::Upper)?;
    //     assert_eq!(manual, linalg);
    //     Ok(())
    // }

    #[test]
    fn quadprogpp_demo() -> Result<()> {
        #[rustfmt::skip]
        let g = array![
            [4.0, -2.0],
            [-2.0, 4.0],
        ];
        let g0 = array![6.0, 0.0];
        #[rustfmt::skip]
        let ce = array![
            [1.0],
            [1.0],
        ];
        let ce0 = array![-3.0];
        #[rustfmt::skip]
        let ci = array![
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        ];
        let ci0 = array![0.0, -2.0, 0.0];
        let mut x = Vector::zeros(g.ncols());
        let best = solve_quadprog(
            g,
            g0,
            ce,
            ce0,
            ci,
            ci0,
            &mut x
        )?;
        assert_ulps_eq!(best, 12.0);
        assert_ulps_eq!(x[0], 1.0);
        assert_ulps_eq!(x[1], 2.0);
        Ok(())
    }

    #[test]
    fn eiquadprog_demo() {
        #[rustfmt::skip]
        let g = array![
            [2.1, 0.0, 1.0],
            [1.5, 2.2, 0.0],
            [1.2, 1.3, 3.1],
        ];
        let g0 = array![6.0, 1.0, 1.0];
        let ce = array![[1.0], [2.0], [-1.0]];
        let ce0 = array![-4.0];
        #[rustfmt::skip]
        let ci = array![
            [1.0, 0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0, -1.0],
            [0.0, 0.0, 1.0,  0.0]
        ];
        let ci0 = array![0.0, 0.0, 0.0, 10.0];
        let mut x = Vector::zeros(g.ncols());
        let best = solve_quadprog(
            g,
            g0,
            ce,
            ce0,
            ci,
            ci0,
            &mut x
        ).unwrap();
        assert_ulps_eq!(best, 6.4);
        assert_ulps_eq!(x[0], 0.0);
        assert_ulps_eq!(x[1], 2.0);
        assert_ulps_eq!(x[2], 0.0);
    }

    // Problem 0 from hmatrix-quadpropp
    #[test]
    fn hmatrix_quadprogpp_problem0() -> Result<()> {
        #[rustfmt::skip]
        let g = array![
            [4.0, 0.0],
            [0.0, 2.0],
        ];
        let g0 = array![-4.0, -8.0];
        #[rustfmt::skip]
        let ci = array![
            [1.0, 0.0, -1.0],
            [0.0, 1.0, -2.0],
        ];
        let ci0 = array![0.0, 0.0, 2.0];
        let ce0 = Vector::zeros(0);
        let ce = Matrix::zeros((g.nrows(), 0));
        let mut x = Vector::zeros(g.ncols());
        let _best = solve_quadprog(
            g,
            g0,
            ce,
            ce0,
            ci,
            ci0,
            &mut x
        )?;
        let answer = x;
        assert_abs_diff_eq!(answer[0], 2.0 / 9.0, epsilon = 1e-12);
        assert_abs_diff_eq!(answer[1], 8.0 / 9.0, epsilon = 1e-12);
        Ok(())
    }

    // Problem 1 from hmatrix-quadprogpp
    #[test]
    fn hmatrix_quadprogpp_problem1() -> Result<()> {
        let offset: Array2<f64> = Array::eye(3) * 1e-12;
        #[rustfmt::skip]
        let g = array![
            [      1.0, 2.0 / 3.0, 1.0 / 3.0],
            [2.0 / 3.0, 2.0 / 3.0,       0.0],
            [1.0 / 3.0,       0.0, 1.0 / 3.0],
        ] + offset;
        let g0 = array![-2.0, -4.0, 2.0];
        let ce = array![[-3.0], [2.0], [1.0]];
        let ce0 = array![0.0];
        #[rustfmt::skip]
        let ci = array![
            [1.0,        0.0,        0.0],
            [0.0,  1.0 / 3.0, -4.0 / 3.0],
            [0.0, -1.0 / 3.0,  1.0 / 3.0]
        ];
        let ci0 = array![0.0, 0.0, 2.0];
        let mut x = Vector::zeros(g.ncols());
        let _best = solve_quadprog(
            g,
            g0,
            ce,
            ce0,
            ci,
            ci0,
            &mut x
        )?;
        let answer = x;
        assert_ulps_eq!(answer[0], 2.0 / 9.0, epsilon = 1e-5);
        assert_abs_diff_eq!(answer[1], 10.0 / 9.0, epsilon = 1e-5);
        assert_abs_diff_eq!(answer[2], -14.0 / 9.0, epsilon = 1e-5);
        Ok(())
    }

}
