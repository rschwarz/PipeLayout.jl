# curried methods for the right-hand side of facts

# can't extend Base.isa, because not generic?!
is_instance(T) = x -> isa(x, T)
