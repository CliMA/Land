macro aimport(mod)
    if isdefined(mod)
        quote
            import $mod
        end
    else
        esc(quote
            arequire(string($(QuoteNode(mod))))
            import $mod
        end)
    end
end

macro ausing(mod)
    if isdefined(mod)
        quote
            using $mod
        end
    else
        esc(quote
            arequire(string($(QuoteNode(mod))))
            using $mod
        end)
    end
end
