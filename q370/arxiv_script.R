library(MASS)
library(tidyverse)
library(jsonlite)
library(igraph)
library(stringr)
library(scales)
library(gtools)

# Colorblind-friendly palette from https://jfly.uni-koeln.de/color/
cb_palette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_fields = c("q-fin" = "#E69F00", "stat" = "#56B4E9", "q-bio" = "#009E73", "eess" = "#F0E442", "cs" = "#0072B2", "physics" = "#D55E00", "math" = "#CC79A7", "total" = "#999999")

## Step 1: Import & clean data
arxiv <- stream_in(file("arxiv-metadata-oai-snapshot.json"))

arxiv_auth <- arxiv %>%
  select(authors, categories) %>%
  mutate(authors = str_replace_all(authors, "\\([^()]*\\)", "")) %>% #remove stuff in parentheses
  mutate(authors = str_replace_all(authors, "\\([^()]*\\)", "")) %>% #remove stuff in parentheses again (b/c nested parentheses)
  mutate(authors = str_replace_all(authors, "\\[|\\]", "")) %>% #remove brackets
  mutate(authors = str_replace_all(authors, "[^:]+: ", "")) %>% #remove stuff before colons (usually a collaboration)
  mutate(authors = str_replace_all(authors, "\\\\'|\\n|\\\\`", "")) %>% #remove formatting artifacts
  mutate(authors = str_replace_all(authors, "^ | $", "")) %>% #remove leading or ending spaces
  mutate(categories = str_extract(categories, "[^ ]+")) #just use the primary category

# Quick tangent to look at author list lengths
# Calculate length of each author list
auth_lengths <- str_split(string=arxiv_auth$authors, ", | and | & |, and ") %>%
  lapply(FUN = length)

auth_lengths <- matrix(unlist(auth_lengths), nrow = length(auth_lengths), byrow = TRUE) %>% 
  data.frame() %>%
  add_column(cat = arxiv_auth$categories) %>%
  rename(length = .)

max(auth_lengths$length) #2829

# Separate author list into individual authors
# Max out at 10 authors
arxiv_auth <- arxiv_auth %>% 
  separate(authors, into = paste0("author", 1:2829), sep = ", | and | & |, and ", remove = FALSE, extra = "drop", fill = "right")

## Step 2: Construct edge list
# Long formatting
arxiv_edges <- arxiv_auth %>%
  filter(!(is.na(author2))) %>% #get rid of papers with one author
  pivot_longer(!c(authors, categories), values_to = "each_author", values_drop_na = TRUE) %>%
  select(!name)

# Take this opportunity to clean some more
arxiv_edges <- arxiv_edges %>%
  filter(!(str_detect(each_author, "The |et al|University|Institute|Collaboration|Center|Tech|Group|Convenors|Contributors|Dept.|Department|Academy"))) %>%
  mutate(each_author = str_replace_all(each_author, "^'|'$", "")) %>%
  mutate(each_author = str_replace_all(each_author, "[^A-Za-z. '-]", "")) %>%
  mutate(each_author = str_replace_all(each_author, "(?<= ) +", "")) %>%
  mutate(each_author = str_replace_all(each_author, "^and ", "")) %>%
  #mutate(each_author = str_replace_all(each_author, "LHCb collaboration ", "")) %>%
  mutate(each_author = str_replace_all(each_author, "^\\.", "")) %>%
  mutate(each_author = str_replace_all(each_author, "^ +| +$", "")) %>%
  filter(!(str_detect(each_author, "^-"))) %>%
  filter(nchar(each_author) > 1)

# Now constructing the combinations between authors for each paper
arxiv_edges <- arxiv_edges %>%
  group_by(authors, categories) %>% 
  filter(n() >= 2) %>% 
  do(data.frame(t(combn(.$each_author, 2)), stringsAsFactors = FALSE)) %>% 
  ungroup() %>%
  select(!authors) %>%
  rename(Source = X1, Target = X2) %>%
  filter(Source != Target)

# Useful function to count nodes
count_nodes <- function(edgelist) {
  return(length(unique(c(unique(edgelist$Source), unique(edgelist$Target)))))
}

# 11,154,834 edges
count_nodes(arxiv_edges) # 1,089,102 nodes

# Collapse edges between the same authors under same category
unique_edges <- arxiv_edges %>%
  rowwise() %>%
  mutate(pair = list(str_sort(c(Source, Target)))) %>% #put each author pair in an alphabetized list so the order is always the same
  ungroup() %>% #remove rowwise operator
  group_by(pair, categories) %>%
  summarize(weight = n()) %>% #weight of the edge is number of co-authored papers
  separate(col = pair, into = c("Source", "Target"), sep = '", "') %>% #separate author pairs again
  mutate(Source = str_replace(Source, '^c\\("', "")) %>% #remove artifacts 
  mutate(Target = str_replace(Target, '"\\)$', ""))

# 7,922,966 edges
count_nodes(unique_edges) # 1,089,102 nodes

# Use just the first listed category (which is generally the most relevant)
cat_edges <- unique_edges %>%
  mutate(cat = str_extract(categories, "[^ ]+")) %>%
  select(!categories) %>%
  mutate(cat = case_when( # collapse aliases
    cat == "solv-int" ~ "nlin.SI",
    cat == "patt-sol" ~ "nlin.PS",
    cat == "adap-org" ~ "nlin.AO",
    cat == "q-alg" ~ "math.QA",
    cat == "alg-geom" ~ "math.AG",
    cat == "dg-ga" ~ "math.DG",
    cat == "cmp-lg" ~ "cs.CL",
    cat == "q-fin.EC" ~ "econ.GN",
    cat == "mtrl-th" ~ "cond-mat.mtrl-sci",
    cat == "chem-ph" ~ "physics.chem-ph",
    cat == "supr-con" ~ "cond-mat.supr-con",
    cat == "comp-gas" ~ "nlin.CG",
    cat == "funct-an" ~ "math.FA",
    cat == "atom-ph" ~ "physics.atom-ph",
    cat == "acc-phys" ~ "physics.acc-ph",
    cat == "plasm-ph" ~ "physics.plasm-ph",
    cat == "ao-sci" ~ "physics.ao-ph",
    cat == "bayes-an" ~ "physics.data-an",
    cat == "chao-dyn" ~ "nlin.CD",
    TRUE ~ cat)) %>%
  group_by(Source, Target, cat) %>%
  summarize(weight = sum(weight)) %>%
  ungroup()

# Use only broad fields
field_edges <- cat_edges %>%
  mutate(cat = str_extract(cat, "[^\\.]+")) %>%
  mutate(cat = case_when(
    str_detect(cat, "^hep|^nucl|-ph$") ~ "physics",
    cat == "cond-mat" ~ "physics",
    cat == "nlin" ~ "physics",
    cat == "gr-qc" ~ "physics",
    cat == "econ" ~ "q-fin",
    TRUE ~ cat))

# Get example edgelist with just connections to Rob
rob <- cat_edges %>% filter(Source == "Robert L. Goldstone" | Target == "Robert L. Goldstone")
write_csv(rob, "rob.csv")

## Step 3: Construct node list
# Get nodes from first column (Source)
nodes1 <- field_edges %>%
  group_by(Source, cat) %>%
  summarize(count = sum(weight))

# Get nodes from second column (Target)
nodes2 <- field_edges %>%
  group_by(Target, cat) %>%
  summarize(count = sum(weight))

# Combine nodes from both columns
nodes <- full_join(nodes1, nodes2, by = c("Source" = "Target", "cat")) %>%
  ungroup()%>%
  mutate(count.x = ifelse(is.na(count.x), 0, count.x)) %>%
  mutate(count.y = ifelse(is.na(count.y), 0, count.y)) %>%
  transmute(author = Source, cat = cat, count = count.x + count.y) %>%
  group_by(author, cat) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  group_by(author) %>%
  filter(count == max(count)) %>%
  mutate(cat = factor(cat, levels = c("q-fin", "eess", "stat", "q-bio", "math", "cs", "physics"), ordered = TRUE)) %>%
  filter(cat == min(cat))

## Step 4: Analysis
# Author lengths
auth_lengths <- auth_lengths %>%
  mutate(cat = case_when( # collapse aliases
    cat == "solv-int" ~ "nlin.SI",
    cat == "patt-sol" ~ "nlin.PS",
    cat == "adap-org" ~ "nlin.AO",
    cat == "q-alg" ~ "math.QA",
    cat == "alg-geom" ~ "math.AG",
    cat == "dg-ga" ~ "math.DG",
    cat == "cmp-lg" ~ "cs.CL",
    cat == "q-fin.EC" ~ "econ.GN",
    cat == "mtrl-th" ~ "cond-mat.mtrl-sci",
    cat == "chem-ph" ~ "physics.chem-ph",
    cat == "supr-con" ~ "cond-mat.supr-con",
    cat == "comp-gas" ~ "nlin.CG",
    cat == "funct-an" ~ "math.FA",
    cat == "atom-ph" ~ "physics.atom-ph",
    cat == "acc-phys" ~ "physics.acc-ph",
    cat == "plasm-ph" ~ "physics.plasm-ph",
    cat == "ao-sci" ~ "physics.ao-ph",
    cat == "bayes-an" ~ "physics.data-an",
    cat == "chao-dyn" ~ "nlin.CD",
    TRUE ~ cat)) %>%
  mutate(cat = str_extract(cat, "[^\\.]+")) %>%
  mutate(cat = case_when(
    str_detect(cat, "^hep|^nucl|-ph$") ~ "physics",
    cat == "cond-mat" ~ "physics",
    cat == "nlin" ~ "physics",
    cat == "gr-qc" ~ "physics",
    cat == "econ" ~ "q-fin",
    TRUE ~ cat))

lengths_df <- auth_lengths %>%
  group_by(cat, length) %>%
  summarize(freq = n()) %>% #get frequency of each length value within each category
  ungroup() %>%
  group_by(cat) %>%
  mutate(prop = freq / sum(freq)) #turn frequency into a proportion

fit <- power.law.fit(auth_lengths$length, xmin=3) #alpha = 2.79

# Plot distribution of author list lengths
ggplot(lengths_df, aes(x = length, y = prop, color = cat)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10(label = label_number()) +
  labs(x = "Length of author list", y = "Proportion of papers") +
  scale_color_manual("Category", values = cb_fields) +
  theme(legend.position = c(0.9, 0.7)) +
  geom_function(fun = function(x) x^(-2.79), xlim = 1:100, color = "black")

edges <- field_edges

# Get frequency of collaborations
collaboration_weights <- edges %>%
  group_by(cat, weight) %>%
  summarize(freq = n()) %>%
  ungroup() %>%
  group_by(cat) %>%
  mutate(prop = freq / sum(freq)) #turn frequency into a proportion

collabs_vector <- edges$weight

fit_collabs <- power.law.fit(collabs_vector, xmin=3)

# Plot edge weight distribution
ggplot(collaboration_weights, aes(x = weight, y = prop, col = cat)) +
  geom_point() +
  scale_y_log10(limits = c(1e-07, 1e+01)) +
  scale_x_log10() +
  labs(x = "Number of co-authored papers", y = "Proportion of collaborations") +
  scale_color_manual("Category", values = cb_fields) +
  theme(legend.position = c(0.9, 0.7)) +
  geom_function(fun = function(x) x^(-2.64), color = "black")

# Get distribution of node degree
g <- graph_from_data_frame(d=edges, directed = FALSE)
d <- degree(g)

fit_degree <- power.law.fit(d, xmin=3) #xmin is important for getting shape right. Only fit the tail of the distribution ignoring first 3

degree_dist <- nodes %>%
  select(cat) %>%
  add_column(degree = d) %>%
  group_by(cat, degree) %>%
  summarize(freq = n()) %>%
  ungroup() %>%
  group_by(cat) %>%
  mutate(prop = freq / sum(freq)) #turn frequency into a proportion

# Plot node degree distribution
ggplot(degree_dist, aes(x = degree, y = prop, color = cat)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10(label = number) +
  labs(y = "Proportion of authors", x = "Number of different co-authors") +
  scale_color_manual("Category", values = cb_fields) +
  theme(legend.position = c(0.9, 0.7)) +
  geom_function(fun = function(x) x^(-1.78), color = "black")

# Edges in each field
edge_field_freq <- edges %>%
  group_by(cat) %>%
  summarize(edges = sum(weight), collabs = n())

# Nodes in each field
#node_field_freq <- nodes %>%
#  group_by(cat) %>%
#  summarize(nodes = n())

node_field_freq <- edges %>%
  select(!weight) %>%
  pivot_longer(cols = !cat, values_to = "authors") %>%
  select(!name) %>%
  group_by(cat) %>%
  summarize(nodes = n_distinct(authors))

mean_total <- lengths_df %>%
  mutate(summed = length * freq) %>%
  ungroup() %>%
  summarize(mean_length = sum(summed)/sum(freq))

mean_auth_lengths <- lengths_df %>%
  mutate(summed = length * freq) %>%
  group_by(cat) %>%
  summarize(mean_length = sum(summed)/sum(freq)) %>%
  add_row(cat = "total", mean_length = mean_total$mean_length)

mean_total <- degree_dist %>%
  mutate(summed = degree * freq) %>%
  ungroup() %>%
  summarize(mean_degree = sum(summed)/sum(freq))

mean_degree <- degree_dist %>%
  mutate(summed = degree * freq) %>%
  group_by(cat) %>%
  summarize(mean_degree = sum(summed)/sum(freq)) %>%
  add_row(cat = "total", mean_degree = mean_total$mean_degree)

mean_total <- collaboration_weights %>%
  mutate(summed = weight * freq) %>%
  ungroup() %>%
  summarize(mean_weight = sum(summed)/sum(freq))

mean_collab_strength <- collaboration_weights %>%
  mutate(summed = weight * freq) %>%
  group_by(cat) %>%
  summarize(mean_weight = sum(summed)/sum(freq)) %>%
  add_row(cat = "total", mean_weight = mean_total$mean_weight)


# Combine into a single data.frame
by_field <- left_join(edge_field_freq, node_field_freq, by = "cat") %>%
  add_row(cat = "total", nodes = sum(.$nodes), edges = sum(.$edges), collabs = sum(.$collabs)) %>% #also throw in the total
  mutate(cat_names = case_when( #make nice names for plotting
    cat == "total" ~ "Total",
    cat == "cs" ~ "Computer Science",
    cat == "math" ~ "Math",
    cat == "stat" ~ "Statistics",
    cat == "q-fin" ~ "Quantitative Finance",
    cat == "q-bio" ~ "Quantitative Biology",
    cat == "eess" ~ "Electrical Engineering\n& Systems Science",
    cat == "physics" ~ "Physics",
  )) #%>%
#arrange(collabs) %>%
#mutate(cat_names = factor(cat_names, levels = .$cat_names, ordered = TRUE)) #so they're plotted in the order of collab frequency

# Get clustering coefficients and giant component size for each field (and total)
g <- graph_from_data_frame(d=edges, directed = FALSE)
component_df <- data.frame(cat = "total", giant_size = max(clusters(g)$csize))
clustering_df <- data.frame(cat = "total", cluster_coef = transitivity(g))
for (i in c("stat", "math", "cs", "q-bio", "q-fin", "eess", "physics")) {
  component_df <- component_df %>% add_row(cat = i, giant_size = max(clusters(graph_from_data_frame(d = filter(edges, cat == i), directed = FALSE))$csize))
  clustering_df <- clustering_df %>% add_row(cat = i, cluster_coef = transitivity(graph_from_data_frame(d = filter(edges, cat == i), directed = FALSE)))
}

# Combine with other data.frame
by_field <- by_field %>%
  left_join(clustering_df, by = "cat") %>%
  left_join(component_df, by = "cat") %>%
  left_join(mean_auth_lengths, by = "cat") %>%
  left_join(mean_degree, by = "cat") %>%
  mutate(giant_percent = giant_size / nodes) %>%
  left_join(mean_collab_strength, by = "cat")

temp <- by_field %>%
  arrange(giant_size) %>%
  mutate(cat_names = factor(cat_names, levels = .$cat_names, ordered = TRUE))

# Plot edges distribution
# excluding total
ggplot(filter(by_field, cat != "total"), aes(x = cat_names, y = edges, color = cat, fill = cat)) +
  #geom_point(size = 5) +
  #geom_segment(aes(xend = cat_names, yend = 0), size = 2) +
  geom_col() +
  #scale_y_log10(labels = label_comma()) +
  scale_y_continuous(labels = label_number_si()) +
  labs(y = "Number of edges", x = "") +
  guides(color = FALSE, fill = FALSE) +
  coord_flip() +
  scale_color_manual(values = cb_fields) +
  scale_fill_manual(values = cb_fields) +
  geom_text(aes(label = comma(edges), y = ifelse(edges > 5000000, edges - 400000, edges + 400000)), color = "black")

# Plot collaborations
ggplot(filter(by_field, cat != "total"), aes(x = cat_names, y = collabs, color = cat, fill = cat)) +
  geom_col() +
  scale_y_continuous(labels = label_number_si()) +
  labs(y = "Number of unique edges (collaborations)", x = "") +
  guides(color = FALSE, fill = FALSE) +
  coord_flip() +
  scale_color_manual(values = cb_fields) +
  scale_fill_manual(values = cb_fields) +
  geom_text(aes(label = comma(collabs), y = ifelse(collabs > 5000000, collabs - 400000, collabs + 400000)), color = "black")

# Plot nodes distribution
ggplot(filter(by_field, cat != "total"), aes(x = cat_names, y = nodes, color = cat, fill = cat)) +
  #geom_point(size = 5) +
  #geom_segment(aes(xend = cat_names, yend = 0), size = 2) +
  geom_col() +
  #scale_y_log10(labels = label_comma()) +
  scale_y_continuous(labels = label_number_si()) +
  labs(y = "Number of nodes", x = "") +
  guides(color = FALSE, fill = FALSE) +
  coord_flip() +
  scale_color_manual(values = cb_fields) +
  scale_fill_manual(values = cb_fields) +
  geom_text(aes(label = comma(nodes), y = ifelse(nodes > 500000, nodes - 40000, nodes + 40000)), color = "black")

# Plot clustering coefficients
ggplot(temp, aes(y = cat_names, x = cluster_coef)) +
  geom_col(aes(fill = cat)) +
  guides(fill = FALSE) +
  scale_fill_manual(values = cb_fields) +
  labs(y = "", x = "Clustering coefficient") +
  geom_text(aes(label = round(cluster_coef, digits = 3), x = ifelse(cluster_coef > 0.55, cluster_coef - 0.03, cluster_coef + 0.03)))

# Plot size of giant component
ggplot(temp, aes(y = cat_names, x = giant_size)) +
  geom_col(aes(fill = cat)) +
  guides(fill = FALSE) +
  scale_fill_manual(values = cb_fields) +
  scale_x_continuous(labels = label_number_si()) +
  labs(y = "", x = "Size of giant component") +
  geom_text(aes(label = comma(giant_size), x = ifelse(giant_size > 750000, giant_size - 75000, giant_size + 70000)))

# Plot size of giant component as percentage
ggplot(temp, aes(y = cat_names, x = giant_percent)) +
  geom_col(aes(fill = cat)) +
  guides(fill = FALSE) +
  scale_fill_manual(values = cb_fields) +
  scale_x_continuous() +
  labs(y = "", x = "Percent of nodes in giant component") +
  geom_text(aes(label = round(giant_percent, 3), x = giant_percent + 0.05))

ggplot(temp, aes(y = cat_names)) +
  geom_col(aes(fill = cat, x = giant_size / giant_percent), alpha = 0.5) +
  geom_col(aes(fill = cat, x = giant_size), alpha = 0.5) +
  guides(fill = FALSE) +
  scale_fill_manual(values = cb_fields) +
  scale_x_continuous(labels = label_number()) +
  labs(y = "", x = "Number of nodes") #+
  #geom_text(aes(label = round(giant_percent, 3), x = giant_size / giant_percent + 0.05))

#g <- graph_from_data_frame(d=edges, directed = FALSE)
#close <- closeness(g)