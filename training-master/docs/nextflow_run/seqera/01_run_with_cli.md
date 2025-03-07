## 1. Use Seqera Platform to capture and monitor Nextflow jobs launched from the CLI

We'll start by using the Nextflow CLI to launch a pipeline and monitor it in Seqera Platform.
Start by logging into the [Seqera Platform](https://cloud.seqera.io/).

!!! info "Nextflow Tower"

    Seqera Platform was previously known as Nextflow Tower.
    You'll still see references to the previous name in environment variables and CLI option names.

### 1.1. Set up your Seqera Platform token by exporting it to your environment

Follow these steps to set up your token:

1.  Create a new token by clicking on the **Settings** drop-down menu:

    ![Create a token](seqera/img/usage_create_token.png)

2.  Name your token:

    ![Name your token](seqera/img/usage_name_token.png)

3.  Save your token safely:

    ![Save token](seqera/img/usage_token.png)

4.  To make your token available to the Nextflow CLI, export it on the command line:

    Open a terminal and type:

    ```bash
     export TOWER_ACCESS_TOKEN=eyxxxxxxxxxxxxxxxQ1ZTE=
    ```

    Where `eyxxxxxxxxxxxxxxxQ1ZTE=` is the token you have just created.

    !!! Warning "Security Note"

        Keep your token secure and do not share it with others.
        Add a ++space++ before the `export` command above to prevent your token from being saved in your shell history.

### 1.2. Run Nextflow CLI with Seqera Platform visualizing and capturing logs

Run a Nextflow workflow with the addition of the `-with-tower` command:

```bash
nextflow run nextflow-io/hello -with-tower
```

You will see output similar to the following:

```console title="Output"
 N E X T F L O W   ~  version 24.04.4

Launching `https://github.com/nextflow-io/hello` [evil_engelbart] DSL2 - revision: afff16a9b4 [master]

Downloading plugin nf-tower@1.9.1
Monitor the execution with Seqera Platform using this URL: https://cloud.seqera.io/user/kenbrewer/watch/5Gs0qqV9Y9rguE
executor >  local (4)
[80/810411] process > sayHello (1) [100%] 4 of 4 ✔
Ciao world!

Bonjour world!

Hola world!

Hello world!
```

Hold ++ctrl++ or ++cmd++ and click on the link to open it in your browser.
You'll see the Seqera Platform interface with the job finished and the logs captured.

![Seqera Platform](seqera/img/run_with_tower.png)

You will see and be able to monitor your **Nextflow jobs** in Seqera Platform.

### 1.3. Set up Seqera Platform in Nextflow configuration

Doing that token setup regularly can become tedious, so let's set this configuration for all our pipeline runs with the global Nextflow configuration file located at `$HOME/.nextflow/config`.

Before we set the configuration, we need to permanently store the token from our environment in Nextflow using a [Nextflow secret](https://www.nextflow.io/docs/latest/secrets.html):

```bash
nextflow secrets set tower_access_token "$TOWER_ACCESS_TOKEN"
```

Make sure your token was saved using:

```bash
nextflow secrets get tower_access_token
```

Next, open the Nextflow configuration file located at `$HOME/.nextflow/config`:

```bash
code $HOME/.nextflow/config
```

Then add the following block of configuration:

```groovy title="$HOME/.nextflow/config"
tower {
    enabled = true
    endpoint = "https://api.cloud.seqera.io"
    accessToken = secrets.tower_access_token
}
```

!!! hint "Endpoint"

    The `endpoint` is the URL of the Seqera Platform API.
    If your institution is running a private instance of Seqera Platform, you should change this to the appropriate URL.

Run your Nextflow workflows as before, but without the `-with-tower` command:

```bash
nextflow run nextflow-io/hello
```

You will see the following output:

```console title="Output"
 N E X T F L O W   ~  version 24.04.4

Launching `https://github.com/nextflow-io/hello` [fabulous_euclid] DSL2 - revision: afff16a9b4 [master]

Monitor the execution with Seqera Platform using this URL: https://cloud.seqera.io/user/kenbrewer/watch/KYjRktIlOuxrh
executor >  local (4)
[71/eaa915] process > sayHello (3) [100%] 4 of 4 ✔
Ciao world!

Bonjour world!

Hola world!

Hello world!
```

Note that we are logging to Seqera Platform even though we did not use the `-with-tower` command!

### 1.4. Use Seqera Platform to explore the resolved configuration of a Nextflow pipeline

Click on the link provided in the output to open the Seqera Platform for your run, then click on the `Configuration` tab.
If you ran your pipeline from the `hello_nextflow` directory, you'll see something like this:

![Seqera Platform Configuration](seqera/img/resolved_configuration.png)

Notice that configuration for our pipeline run is being run pulled from three separate files:

- `/home/gitpod/.nextflow/config` - This is the global configuration file we just added.
- `/home/gitpod/.nextflow/assets/nextflow-io/hello/nextflow.config` - This is the `nextflow.config` file from the `nextflow-io/hello` repository.
- `/workspaces/training/nf-training/hello-nextflow/nextflow.config` - This is the `nextflow.config` file from our current working directory.

Nextflow resolves these configurations at runtime with a [specific order of precedence](https://www.nextflow.io/docs/latest/config.html#configuration-file).
The general rule, however, is that more specific configurations override less specific ones, and config/params specified on the CLI will override defaults in the config files.

Helpfully, Seqera Platform shows us the final output of this configuration resolution process which can be very useful for debugging!

### Takeaway

You have learned how to:

- Set up your Seqera Platform token by exporting it to your environment.
- Run Nextflow CLI with Seqera Platform visualizing and capturing logs.
- Set up Seqera Platform logging by default.
- Use Seqera Platform to explore the resolved configuration of a Nextflow pipeline.

### What's next?

Learn how to launch Nextflow pipelines from Seqera Platform using the Launchpad feature.
