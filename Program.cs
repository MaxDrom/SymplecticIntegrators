using System.Globalization;
using System.Numerics;

namespace SymplecticIntegrators;

class Program
{

    static decimal Newton(Vector<decimal> q)
    {
        var r = MyMath.Root(q * q, 2);
        return 1 / r;
    }

    static Vector<decimal> dNewton(Vector<decimal> q)
    {
        var r = MyMath.Root(q * q, 2);
        r = 1 / (r * r * r);
        return q * r;
    }

    static void Main(string[] args)
    {
        CultureInfo.CurrentCulture = new CultureInfo("en-US", false);
        CultureInfo.CurrentCulture.NumberFormat.CurrencyDecimalDigits = 28;

        var taus = new decimal[] { 0.1M, 0.01M, 0.001M, 0.0001M };
        var T = 10M;
        var q0 = new Vector<decimal>([1.0000000000000000000000000000M]);
        var p0 = new Vector<decimal>([0.0000000000000000000000000000M]);
        foreach (var tau in taus)
            CreateData(tau, T, z => z, z => z, (q, p) => (q * q + p * p) / 2, q0, p0, "pendulum");

        q0 = new Vector<decimal>([1.0000000000000000000000000000M, 0.0000000000000000000000000000M]);
        p0 = new Vector<decimal>([0.0000000000000000000000000000M, 1.0000000000000000000000000000M]);
        foreach (var tau in taus)
            CreateData(tau, T, dNewton, z => z, (q, p) => -Newton(q) + p * p / 2, q0, p0, "orbit");

        var yoshida = YoshidaIntegrator<double, Vector<double>>.BuildFromLeapfrog(z=>new Vector<double>(z.Select(z=>Math.Sin(z))), z=>z, 2);
        foreach(var x0 in new double[]{-6, -4, -3, -2, -1} )
        {
            using var file = new StreamWriter($"results/{x0}_yoshida_idpen.dat");
            
            foreach (var r in yoshida.Integrate(100.0, 0.01, new Vector<double>([x0]), new Vector<double>([0.0])))
            {
                
                var (t, q, p) = r;
                var normalizedDeg = q[0]%(2*Math.PI);
                if (normalizedDeg <= -Math.PI)
                    normalizedDeg += 2*Math.PI;
                else if (normalizedDeg > Math.PI)
                    normalizedDeg -= 2*Math.PI;
                file.WriteLine(normalizedDeg + " " + p + " ");
            }
        }

            using var file1 = new StreamWriter($"results/{-2.5}_yoshida_idpen.dat");
            
            foreach (var r in yoshida.Integrate(100.0, 0.01, new Vector<double>([-2.5]), new Vector<double>([2.0])))
            {
                
                var (t, q, p) = r;
                var normalizedDeg = q[0]%(2*Math.PI);
                if (normalizedDeg <= -Math.PI)
                    normalizedDeg += 2*Math.PI;
                else if (normalizedDeg > Math.PI)
                    normalizedDeg -= 2*Math.PI;
                file1.WriteLine(normalizedDeg + " " + p + " ");
            }

            using var file2 = new StreamWriter($"results/{2.5}_yoshida_idpen.dat");
            
            foreach (var r in yoshida.Integrate(100.0, 0.01, new Vector<double>([-2.5]), new Vector<double>([-2.0])))
            {
                
                var (t, q, p) = r;
                var normalizedDeg = q[0]%(2*Math.PI);
                if (normalizedDeg <= -Math.PI)
                    normalizedDeg += 2*Math.PI;
                else if (normalizedDeg > Math.PI)
                    normalizedDeg -= 2*Math.PI;
                file2.WriteLine(normalizedDeg + " " + p + " ");
            }
    }
    
    static void CreateData<TField>(TField tau, TField T, Func<Vector<TField>, Vector<TField>> dV,
            Func<Vector<TField>, Vector<TField>> dQ,
            Func<Vector<TField>, Vector<TField>, TField> H,
            Vector<TField> q0, Vector<TField> p0, string systemName)
        where TField: IFloatingPoint<TField>
    {
        var integrators = new Dictionary<string, SymplecticIntegrator<TField, Vector<TField>>>();
        integrators["euler"] = new SymplecticEuler<TField, Vector<TField>>(dV, dQ);
        integrators["leapfrog"] = new Leapfrog<TField, Vector<TField>>(dV, dQ);
        integrators["yoshida"] = YoshidaIntegrator<TField, Vector<TField>>.BuildFromLeapfrog(dV, dQ, 2);
        integrators["yoshida6"] = YoshidaIntegrator<TField, Vector<TField>>.BuildFromLeapfrog(dV, dQ, 3);
        foreach (var integrator in integrators)
        {
            using var file = new StreamWriter($"results/{tau}_{integrator.Key}_{systemName}.dat");
            
            foreach (var r in integrator.Value.Integrate(T, tau, q0, p0))
            {
                
                var (t, q, p) = r;
                var e = H(q, p);
                
                file.WriteLine(t + " " + q + " " + p + " " + e);
            }
        }
    }
}

