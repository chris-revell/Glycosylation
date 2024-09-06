
module CisternaWidth

using XLSX
using DataFrames
using DrWatson
using Interpolations

function hFromData(dimsPlus; cisternaSeriesID=1)
    df = DataFrame(XLSX.readtable(datadir("exp_raw", "ModellingData.xlsx"), 1))
    hs = filter(x->x.series_ID==cisternaSeriesID, df)[1,5:end]
    hs = collect(skipmissing(hs))./100.0
    xMaxInferred = (length(hs)-1)   
    xs   = collect(range(0.0, xMaxInferred, dimsPlus[2])) # Positions of discretised vertices in space
    itp_cubic = cubic_spline_interpolation(0:xMaxInferred, hs)
    hFun(x) = itp_cubic(x)
    mat_h = zeros(dimsPlus)
    for j=1:dimsPlus[1]
        # mat_h[:, j] .= hFun(xs[j])
        selectdim(mat_h, 2, j) .= hFun(xs[j])
    end
    return xMaxInferred, mat_h, xs
end

function hFromFunction(dimsPlus; xMax=100.0, μxh = 50.0, σxh = 10.0)
    hFun(x) = 0.1+exp(-(x-μxh)^2/σxh^2)
    xs   = collect(range(0.0, xMax, dimsPlus[2])) # Positions of discretised vertices in space
    mat_h = zeros(dimsPlus...)
    for j=1:dimsPlus[2]
        # mat_h[:, j] .= hFun(xs[j])
        selectdim(mat_h, 2, j) .= hFun(xs[j])
    end
    return xMax, mat_h, xs
end

export hFromData, hFromFunction

end

# function hFromFunction(dimsPlus; xMax=100.0, μxh = 50.0, σxh = 10.0)
#     hFun(x) = 0.1+exp(-(x-μxh)^2/σxh^2)
#     xs   = collect(range(0.0, xMax, dimsPlus[2])) # Positions of discretised vertices in space
#     mat_h = zeros(dimsPlus)
#     if length(dimsPlus)==2
#         for jj=1:dimsPlus[2]
#             # mat_h[:, j] .= hFun(xs[j])
#             selectdim(mat_h, 2, jj) .= 0.1+exp(-(xs[jj]-μxh)^2/σxh^2) #hFun(xs[jj])
#         end
#     else
#         for kk=1:dimsPlus[3], jj=1:dimsPlus[2]
#             mat_h[:, jj, kk] .= 0.1+exp(-(xs[jj]-μxh)^2/σxh^2 - )
#             # selectdim(mat_h, 2, jj) .= hFun(xs[jj])
#         end
#     end
#     return xMax, mat_h, xs
# end

# function hFromFunction(dimsPlus; xMax=100.0, μxh = 50.0, σxh = 10.0)
#     hFun(x) = 0.1+exp(-(x-μxh)^2/σxh^2)
#     xs   = collect(range(0.0, xMax, dimsPlus[2])) # Positions of discretised vertices in space
#     mat_h = zeros(dimsPlus)
#     if length(dimsPlus)==2
#         for jj=1:dimsPlus[2]
#             # mat_h[:, j] .= hFun(xs[j])
#             selectdim(mat_h, 2, jj) .= 0.1+exp(-(xs[jj]-μxh)^2/σxh^2) #hFun(xs[jj])
#         end
#     else
#         for kk=1:dimsPlus[3], jj=1:dimsPlus[2]
#             mat_h[:, jj, kk] .= 0.1+exp(-(xs[jj]-μxh)^2/σxh^2 - )
#             # selectdim(mat_h, 2, jj) .= hFun(xs[jj])
#         end
#     end
#     return xMax, mat_h, xs
# end