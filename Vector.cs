using System.Collections;
using System.Numerics;

namespace SymplecticIntegrators;

public class Vector<TField> : ILinearSpace<Vector<TField>, TField>, IEnumerable<TField>
    where TField : INumber<TField>

{
    private bool _isZero = false;
    private TField[] _coords;

    private Vector(bool isZero = false)
    {
        _coords = [];
        _isZero = isZero;
    }
    public Vector(int n)
    {
        _coords = new TField[n];
    }

    public Vector(IEnumerable<TField> coords)
    {
        _coords = coords.ToArray();
    }

    public TField this[int index]
    {
        get => _coords[index];
        set => _coords[index] = value;
    }

    public static Vector<TField> AdditiveIdentity => new Vector<TField>(true);

    public int Length => _coords.Length;

    public static Vector<TField> operator +(Vector<TField> left, Vector<TField> right)
    {
        if (left._isZero)
            return right;
        if (right._isZero)
            return left;
#if DEBUG
        if (left.Length != right.Length)
            throw new ArgumentException("Вектора должны иметь одинаковый размер");
#endif
        var result = new Vector<TField>(left.Length);
        for (var i = 0; i < left.Length; i++)
            result[i] = left[i] + right[i];
        return result;
    }

    public static Vector<TField> operator *(Vector<TField> left, TField right)
    {
        if (left._isZero)
            return new Vector<TField>(true);
        var result = new Vector<TField>(left.Length);
        for (var i = 0; i < left.Length; i++)
            result[i] = left[i] * right;
        return result;
    }

    public static Vector<TField> operator *(TField right, Vector<TField> left)
    {
        return left * right;
    }

    public static TField operator *(Vector<TField> right, Vector<TField> left)
    {
        var result = TField.Zero;
        if(left._isZero || right._isZero)
            return TField.Zero;
#if DEBUG
        if (left.Length != right.Length)
            throw new ArgumentException("Вектора должны иметь одинаковый размер");
#endif
        for(var i = 0; i<left.Length; i++)
            result += right[i]*left[i];

        return result;
    }

    public static Vector<TField> operator -(Vector<TField> value)
    {
        var result = new Vector<TField>(value.Length);
        for(var i = 0; i<value.Length; i++)
            result[i] = -value[i];

        return result;
    }

    public static Vector<TField> operator -(Vector<TField> left, Vector<TField> right)
    {
        return left +(-right);
    }

    public override string ToString()
    {
        return string.Join(" ", _coords);
    }

    public IEnumerator<TField> GetEnumerator()
    {
        return ((IEnumerable<TField>)_coords).GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return _coords.GetEnumerator();
    }
}