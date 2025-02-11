"""
Custom type: year
"""
struct Year
    value::Int
end
Base.:(==)(a::Year, b) = Base.:(==)(a.value, b)
Base.:(==)(a, b::Year) = Base.:(==)(a, b.value)
Base.:(==)(a::Year, b::Year) = Base.:(==)(a.value, b.value)
Base.:(+)(a::Year, b) = Year(Base.:(+)(a.value, b))
Base.show(io::IO, data::Year) = Base.show(io::IO, data.value)
Base.isless(a::Year, b::Year) = Base.isless(a.value, b.value)
Base.isless(a::Year, b) = Base.isless(a.value, b)
Base.isless(a, b::Year) = Base.isless(a, b.value)

"""
Custom type: peer
"""
struct Peer
    value::String
end
Base.:(==)(a::Peer, b) = Base.:(==)(a.value, b)
Base.:(==)(a, b::Peer) = Base.:(==)(a, b.value)
Base.:(==)(a::Peer, b::Peer) = Base.:(==)(a.value, b.value)
Base.show(io::IO, data::Peer) = Base.show(io::IO, data.value)
Base.isless(a::Peer, b::Peer) = Base.isless(a.value, b.value)
Base.isless(a::Peer, b) = Base.isless(a.value, b)
Base.isless(a, b::Peer) = Base.isless(a, b.value)

"""
Custom type: technology
"""
struct Tec
    value::String
end
Base.:(==)(a::Tec, b) = Base.:(==)(a.value, b)
Base.:(==)(a, b::Tec) = Base.:(==)(a, b.value)
Base.:(==)(a::Tec, b::Tec) = Base.:(==)(a.value, b.value)
Base.show(io::IO, data::Tec) = Base.show(io::IO, data.value)
Base.isless(a::Tec, b::Tec) = Base.isless(a.value, b.value)
Base.isless(a::Tec, b) = Base.isless(a.value, b)
Base.isless(a, b::Tec) = Base.isless(a, b.value)
