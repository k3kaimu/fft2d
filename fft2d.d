// Written in the D programming language

/**
2D FFT Library

License:    NYSL

Authors:    Kazuki Komatsu
*/

module fft2d;

import std.algorithm,
       std.complex,
       std.conv,
       std.exception,
       std.range,
       std.traits,
       std.typetuple;


template isComplexLike(T)
{
    enum isComplexLike = is(typeof((T t){
            static assert(isFloatingPoint!(typeof(t.re)));
            static assert(isFloatingPoint!(typeof(t.im)));
        }));
}


template isRangeOfRanges(alias pred, T)
{
    enum isRangeOfRanges = pred!T && pred!(ElementType!T);
}


void main()
{
    import std.algorithm, std.random, std.range, std.numeric, std.math;

    auto obj = new Fft2D!Fft(16, 1);
    auto input = [iota(0, 16).map!(a => sin((a*2/16.0)*PI))];
    auto result = obj.fft!real(input);
    auto inv = obj.inverseFft!real(result);

    import std.stdio;
    writeln(input);
    writeln(result);
    writeln(inv);
}


class Fft2D(FFTObj)
{
    this(int x, int y)
    {
        _xfftObj = new FFTObj(x);
        _yfftObj = new FFTObj(y);
    }


    @property size_t xsize() const
    {
        return _xfftObj.size;
    }


    @property size_t ysize() const
    {
        return _yfftObj.size;
    }


    Complex!F[][] fft(F = double, RofR)(RofR range) const
    {
        return fftImpl!("fft", F)(range);
    }


    void fft(Ret, RofR)(RofR range, Ret buf) const
    {
        fftImpl!"fft"(range, buf);
    }


    Complex!F[][] inverseFft(F = double, RofR)(RofR range) const
    {
        return fftImpl!("inverseFft", F)(range);
    }


    void inverseFft(Ret, RofR)(RofR range, Ret buf) const
    {
        fftImpl!"inverseFft"(range, buf);
    }



  private:
    void enforceSize(RofR)(RofR range) const
    {
        enforce(range.length <= ysize, text("FFT size mismatch.  Expected ", [ysize, xsize], ", got ", [range.length, range[0].length]));

        foreach(e; range)
            enforce(e.length <= xsize, text("FFT size mismatch.  Expected ", [ysize, xsize], ", got ", [range.length, e.length]));
    }


    Complex!F[][] fftImpl(string op = "fft", F, RofR)(RofR range) const
    if(isFloatingPoint!F
    && isRangeOfRanges!(isRandomAccessRange, RofR) && (isFloatingPoint!(ElementType!(ElementType!RofR)) || isComplexLike!(ElementType!(ElementType!RofR))))
    {
        enforceSize(range);
        Complex!F[][] ret;
        
        if(range.length == 0 || range[0].length == 0)
            return ret;

        ret = uninitializedArray!(Complex!F[][])(ysize, xsize);
        
        mixin("this." ~ op ~ "(range, ret);");
        return ret;
    }


    void fftImpl(string op = "fft", Ret, RofR)(RofR range, Ret buf) const
    if(isRangeOfRanges!(templateAnd!(isRandomAccessRange, hasSlicing), RofR) && (isFloatingPoint!(ElementType!(ElementType!RofR)) || isComplexLike!(ElementType!(ElementType!RofR)))
    && isRangeOfRanges!(templateAnd!(isRandomAccessRange, hasSlicing), Ret) && hasLvalueElements!(ElementType!Ret) && isComplexLike!(ElementType!(ElementType!Ret)))
    {
        alias F = typeof(ElementType!(ElementType!Ret).init.re);

        enforceSize(range);

        immutable xsize = this.xsize,
                  ysize = this.ysize,
                  bufSizeByte = xsize * ysize * Complex!F.sizeof;

        if(xsize == 0 || ysize == 0)
            return;

        // resize _buffer
        if(_buffer.length < bufSizeByte)
            _buffer.length = bufSizeByte;

        scope cpxBuffer = () @trusted { return cast(Complex!F[])_buffer; }();

        foreach(i; 0 .. ysize)
            mixin("_xfftObj." ~ op ~ "(range[i], cpxBuffer[i*xsize .. (i+1)*xsize]);");


        static struct ByIndex2D()
        {
            @property auto ref front() { return _r.front[_idx]; }
            @property auto ref back() { return _r.back[_idx]; }
            @property bool empty() { return _r.empty; }
            void popFront() { _r.popFront(); }
            void popBack() { _r.popBack(); }
            @property typeof(this) save() { return typeof(this)(_r.save, _idx); }
            typeof(this) opSlice(size_t i, size_t j) { return typeof(this)(_r[i .. j], _idx); }
            @property size_t length() { return _r.length; }
            alias length opDollar;
            auto ref opIndex(size_t i) { return _r[i][_idx]; }

          private:
            Ret _r;
            size_t _idx;
        }

        // y方向、つまり各列でのフーリエ変換
        foreach(i; 0 .. xsize)
            mixin("_yfftObj." ~ op ~ "(cpxBuffer[i .. $].stride(xsize).take(ysize), ByIndex2D!()(buf, i));");
    }


    FFTObj _xfftObj;
    FFTObj _yfftObj;
    
  static:
    void[] _buffer;
}
