using System.Numerics;

namespace SymplecticIntegrators;

public class YoshidaIntegrator<TField, TSpace> : SymplecticIntegrator<TField, TSpace>
    where TField : IFloatingPoint<TField>
    where TSpace : ILinearSpace<TSpace, TField>
{
    private SymplecticIntegrator<TField, TSpace> _previos;
    private TField _x0;
    private TField _x1;

    public int Order {get; private set;}

    private (TField, TField) FindX(int order)
    {
        Order = order;
        var x1 =  TField.One/(_two - MyMath.Root(_two, 2*order +1));
        var x0 = TField.One - _two*x1;
 
        var fieldorder = (dynamic) order;

        for (var i = 0; i<100; i++)
        {
            var w11 = TField.One;
            var w12 = _two;
            var w21 = _two*MyMath.Pow(x1, 2*order )*(_two*fieldorder + TField.One);
            var w22 = MyMath.Pow(x0, 2*order )*(_two*fieldorder + TField.One);
            var det = w11*w22 - w21*w12;
            var y0 = _two*x1 +x0-TField.One;
            var y1 = _two*MyMath.Pow(x1, 2*order +1) + MyMath.Pow(x0, 2*order +1);

            x0 = x0-(w22*y0-w12*y1)/det;
            x1 = x1-(-w21*y0+w11*y1)/det;
        }
        return (x0, x1);
    }

    private YoshidaIntegrator(int order,
        SymplecticIntegrator<TField, TSpace> previosIntegrator) : base(previosIntegrator.dV, previosIntegrator.dT)
    {
        _previos = previosIntegrator;

        (_x0, _x1) = FindX(order);
    }

    public override (TSpace, TSpace) Step(TSpace q0, TSpace p0, TField tau)
    {
        var (q, p) = _previos.Step(q0, p0, tau * _x1);
        (q, p)= _previos.Step(q, p, tau * _x0);
        return _previos.Step(q, p, tau * _x1);
    }

    public static SymplecticIntegrator<TField, TSpace> BuildFromLeapfrog(Func<TSpace, TSpace> dV, Func<TSpace, TSpace> dT, int order)
    {
        return BuildFromBasic(new Leapfrog<TField, TSpace>(dV, dT), order);
    }

    private static SymplecticIntegrator<TField, TSpace> BuildFromBasic(SymplecticIntegrator<TField, TSpace> basic, int order)
    {
        if (order == 1)
            return basic;
        var previosIntegrator = basic;
        for (int i = 1; i < order; i++)
        {
            previosIntegrator = new YoshidaIntegrator<TField, TSpace>(i, basic);
        }
        return previosIntegrator;
    }
}