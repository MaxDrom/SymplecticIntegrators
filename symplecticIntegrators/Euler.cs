using System.Numerics;

namespace SymplecticIntegrators;

public class SymplecticEuler<TField, TSpace> : SymplecticIntegrator<TField, TSpace>
    where TField : IFloatingPoint<TField>
    where TSpace : ILinearSpace<TSpace, TField>
{
    public SymplecticEuler(Func<TSpace, TSpace> dV, Func<TSpace, TSpace> dT) : base(dV, dT)
    {
    }

    public override (TSpace, TSpace) Step(TSpace q0, TSpace p0, TField tau)
    {
        var result = StepByV(tau, q0, p0);
        result = StepByT(tau, result.Item1, result.Item2);
        return result;
    }
}