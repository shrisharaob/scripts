function y = ConvertFs(x, inFs, outFs)

    y = round(x .* outFs ./ inFs) + 1;

end