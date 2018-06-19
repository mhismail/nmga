project_proto <- list(name = character(),
                      path = character(),
                      control_stream = character(),
                      data = data.frame(),
                      tokens = data.frame(),
                      models = list(),
                      analyses = list()
)

slot.names <- names(project_proto)
slots <- sapply(project_proto, class)
names(slots) <- names(project_proto)

setClass("project", slots = slots, prototype = project_proto)
