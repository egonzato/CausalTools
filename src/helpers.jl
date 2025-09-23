module Helpers
export check_link

using GLM
using Logging   # for @warn

function check_link(link)
    allowed_links = [LogitLink, ProbitLink]
    if !any(isa(link, L) for L in allowed_links)
        @warn "The link function you provided is not allowed. Allowed: $(allowed_links)"
        error("Invalid link function. Execution stopped.")
    end
end

end
