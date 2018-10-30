dedup <- read_tsv("clusters/dedup.tsv") %>%
  rename(dedup = cluster, gene = member)

dedup_count <- dedup %>%
  group_by(dedup) %>%
  summarize(dedup_count = n())

cascade <- read_tsv("clusters/cascade.tsv") %>%
  rename(cascade = cluster, dedup = member)

cascade_count <- cascade %>%
  group_by(cascade) %>%
  summarize(cascade_count = n())

profile <- read_tsv("clusters/profile.tsv") %>%
  rename(profile = cluster, dedup = member)

profile_count <- profile %>%
  group_by(profile) %>%
  summarize(profile_count = n())

joined <- dedup %>%
  left_join(cascade, by="dedup") %>%
  left_join(profile, by="dedup") %>%
  left_join(dedup_count, by="dedup") %>%
  left_join(cascade_count, by="cascade") %>%
  left_join(profile_count, by="profile") %>%
  select(gene, dedup, dedup_count, cascade, cascade_count, profile, profile_count) %>%
  arrange(
    desc(profile_count),
    profile,
    desc(cascade_count),
    cascade,
    desc(dedup_count),
    dedup,
    gene
  )

write_tsv(joined, "joined.tsv")
