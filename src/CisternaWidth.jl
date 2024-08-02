
module CisternaWidth

using XLSX
using DataFrames
using DrWatson
using Interpolations

function hFunFromData(; cisternaSeriesID=1)
    df = DataFrame(XLSX.readtable(datadir("exp_raw", "ModellingData.xlsx"), 1))
    hs = filter(x->x.series_ID==cisternaSeriesID, df)[1,5:end]
    hs = collect(skipmissing(hs))./100.0
    xMaxInferred = (length(hs)-1)   
    itp_cubic = cubic_spline_interpolation(0:xMaxInferred, hs)
    hFunD(x) = itp_cubic(x)
    return xMaxInferred, hFunD
end

function hFunFromFunction(; xMax=100.0)
    μxh = xMax/2.0; σxh = xMax/10.0
    # hFun(x, μx, σx) = 0.1+exp(-(x-μx)^2/σx^2)
    hFunF(x) = 0.1+exp(-(x-μxh)^2/σxh^2)
    return xMax, hFunF
end

export hFunFromData, hFunFromFunction

end
