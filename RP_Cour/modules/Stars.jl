function mystars(x)
    if x < 0.001
        y = "***"
    elseif x < 0.05
        y = "**"
    elseif x < 0.10
        y =  "*"
    else
        y = ""
    end
end

# mystars.(table[:,2])
