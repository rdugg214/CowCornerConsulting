function result = CalculateJacobian(type, i, m, err, preverr)
    if (type == "Full Newton")
        result = true;
    elseif (type == "Chord")
        result = false;
    else
        result = mod(i,m)==0 || preverr < err;
    end
end